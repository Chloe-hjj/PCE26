import numpy as np
import pandas as pd
import freesasa
import yaml
from string import ascii_uppercase
import MDAnalysis as mda
import argparse
import sys

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

amino_acid_properties = {
    "ALA": "Nonpolar",
    "VAL": "Nonpolar",
    "LEU": "Nonpolar",
    "ILE": "Nonpolar",
    "PHE": "Nonpolar",
    "TRP": "Nonpolar",
    "MET": "Nonpolar",
    "PRO": "Nonpolar",
    "GLY": "Nonpolar",
    "SER": "Polar",
    "THR": "Polar",
    "CYS": "Polar",
    "TYR": "Polar",
    "ASN": "Polar",
    "GLN": "Polar",
    "LYS": "Basic",
    "ARG": "Basic",
    "HIS": "Basic",
    "ASP": "Acidic",
    "GLU": "Acidic"
}

# Functions for identifying candidate GLU/ASP and LYS
def compute_sasa(
    pdb_path: str
) -> dict[str, list[dict]]:
    """Return exposed GLU/ASP and LYS residues sorted by SASA (descending)."""
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    residue_areas = result.residueAreas()

    out = []
    
    for chain_label, residues in residue_areas.items():
        for res_number, area in residues.items():
            
            entry = {
                "chain": chain_label,
                "residue_number": res_number,
                "residue_type": area.residueType,
                "sasa_total": round(area.total, 2),
                "sasa_side_chain": round(area.sideChain, 2),
                "relative_side_chain": round(area.relativeSideChain, 3),
            }

            out.append(entry)

    return pd.DataFrame.from_records(out)

def compute_distance_to_center(
    pdb_path: str
) -> dict[str, list[dict]]:
    """Return exposed GLU/ASP and LYS residues sorted by SASA (descending)."""
    structure = mda.Universe(pdb_path)
    positions = structure.atoms.select_atoms("name NZ").positions
    distances = np.linalg.norm(positions, axis=1)
    resnums = structure.atoms.select_atoms("name NZ").resnums
    resnames = structure.atoms.select_atoms("name NZ").resnames
    chids = structure.atoms.select_atoms("name NZ").chainIDs
    u = pd.DataFrame.from_dict(
        dict(chain=chids, residue_number=resnums, residue_type=resnames, distance=distances)
    )
    u.residue_number = u.residue_number.astype(str)
    return u

def extract_chains(pdb_path: str) -> dict:
    u = mda.Universe(pdb_path)
    protein = u.select_atoms("protein")
    result = {}
    for chain_id in dict.fromkeys(protein.chainIDs):
        residues = protein.select_atoms(f"chainID {chain_id}").residues
        seq = "".join(THREE_TO_ONE.get(r.resname, "X") for r in residues)
        res_map = {str(r.resid): i for i, r in enumerate(residues, 1)}
        result[chain_id] = (seq, res_map)
    return result

def main():
    parser = argparse.ArgumentParser(description="Generate Boltz-2 input")
    parser.add_argument("--pdb1",  required=True, help="Path to first PDB (Lysine source)")
    parser.add_argument("--pdb2",  required=True, help="Path to second PDB (Acidic source)")
    parser.add_argument("--output", default="boltz_input.yaml")
    args = parser.parse_args()
    
    pdb1 = args.pdb1
    pdb2 = args.pdb2
    
    # Analysis for pdb1 (Basic source)
    sasa_df1 = compute_sasa(pdb1)
    dist_df1 = compute_distance_to_center(pdb1)

    sasa_df1['residue_family'] = sasa_df1['residue_type'].map(amino_acid_properties)

    # Merge and normalize distances for pdb1
    full_df1 = pd.merge(dist_df1, sasa_df1)

    full_df1['distance_mean'] = full_df1['distance'].mean()
    full_df1['distance_norm'] = full_df1['distance'] - full_df1['distance_mean']
    full_df1['location'] = full_df1['distance_norm'].apply(lambda x: 'outer' if x > 0 else 'inner')

    # Identify the best LYS
    best_lys = (full_df1.query('location == "outer" and residue_type == "LYS"').sort_values('sasa_total', ascending = False).head(1))

    # Analysis for pdb2 (Acidic source) - Now tracking SASA only
    sasa_df2 = compute_sasa(pdb2)
    sasa_df2['residue_family'] = sasa_df2['residue_type'].map(amino_acid_properties)

    # We skip the distance calculations entirely and use the SASA data directly
    full_df2 = sasa_df2

    # Identify the best GLU or ASP based strictly on the highest SASA score
    best_acid = (full_df2.query('residue_type == "GLU" or residue_type == "ASP"')
                 .sort_values('sasa_total', ascending=False)
                 .head(1))

    # Extract info and determine the reactive atom
    acid_type = best_acid['residue_type'].values[0]
    acid_atom = "CD" if acid_type == "GLU" else "CG"

    chains1 = extract_chains(pdb1)
    chains2 = extract_chains(pdb2)

    boltz_id_iter = iter(ascii_uppercase)
    pdb1_mapping = {c: next(boltz_id_iter) for c in chains1}
    pdb2_mapping = {c: next(boltz_id_iter) for c in chains2}

    sequences = []
    for cid, (seq, _) in chains1.items():
        sequences.append({"protein": {"id": pdb1_mapping[cid], "sequence": seq}})
    for cid, (seq, _) in chains2.items():
        sequences.append({"protein": {"id": pdb2_mapping[cid], "sequence": seq}})

    # Map PDB numbers to Boltz Sequence indices
    lys_chain = best_lys['chain'].values[0]
    lys_idx = chains1[lys_chain][1][str(best_lys['residue_number'].values[0])]

    acid_chain = best_acid['chain'].values[0]
    acid_idx = chains2[acid_chain][1][str(best_acid['residue_number'].values[0])]
    
    constraints = [{
        "bond": {
            "atom1": [pdb1_mapping[lys_chain], lys_idx, "NZ"],
            "atom2": [pdb2_mapping[acid_chain], acid_idx, acid_atom],
        }
    }]
   
    boltz_input = {"version": 1, "sequences": sequences, "constraints": constraints}
    with open(args.output, "w") as f:
        yaml.dump(boltz_input, f, default_flow_style=False, sort_keys=False)

    print(f"Successfully generated {args.output}")
    
if __name__ == "__main__":
    main()
