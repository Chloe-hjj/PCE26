#!/usr/bin/env python3
"""Generate a Boltz-2 input YAML with a covalent bond between the top exposed
lysine of protein 1 and the top exposed carboxyl (GLU/ASP) of protein 2.

Reads the YAML produced by sasa_exposed_residues.py and the two PDB files to
extract chain sequences and map PDB residue numbers to 1-based sequence indices
required by Boltz-2.

The generated YAML is intended for use with ``boltz predict --use_msa_server``.
"""

from string import ascii_uppercase

import click
import MDAnalysis as mda
import yaml

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

LYS_ATOM = "NZ"
CARBOXYL_ATOM = {"GLU": "CD", "ASP": "CG"}


def extract_chains(pdb_path: str) -> dict[str, tuple[str, dict[str, int]]]:
    """Return ``{chain: (sequence, {pdb_resnum_str: 1-based_idx})}`` from *pdb_path*."""
    u = mda.Universe(pdb_path)
    protein = u.select_atoms("protein")

    result: dict[str, tuple[str, dict[str, int]]] = {}
    for chain_id in dict.fromkeys(protein.chainIDs):
        residues = protein.select_atoms(f"chainID {chain_id}").residues
        seq = "".join(THREE_TO_ONE.get(r.resname, "X") for r in residues)
        res_map = {str(r.resid): i for i, r in enumerate(residues, 1)}
        result[chain_id] = (seq, res_map)
    return result


@click.command()
@click.argument("sasa_yaml", type=click.Path(exists=True, dir_okay=False))
@click.argument("pdb1", type=click.Path(exists=True, dir_okay=False))
@click.argument("pdb2", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o", "--output", default="boltz_input.yaml", show_default=True,
    help="Output Boltz-2 YAML file.",
)
def main(sasa_yaml: str, pdb1: str, pdb2: str, output: str) -> None:
    """Generate a Boltz-2 input YAML with a covalent bond constraint.

    SASA_YAML is the output of sasa_exposed_residues.py.
    PDB1 / PDB2 are the two protein structures (order must match the SASA input).
    """
    # --- Load SASA results -------------------------------------------------
    with open(sasa_yaml) as fh:
        sasa = yaml.safe_load(fh)

    struct1, struct2 = sasa["structures"][0], sasa["structures"][1]

    if not struct1["basic_exposed"]:
        raise click.ClickException("No exposed lysines found in protein 1.")
    if not struct2["acidic_exposed"]:
        raise click.ClickException("No exposed carboxyls found in protein 2.")

    top_lys = struct1["basic_exposed"][0]
    top_acid = struct2["acidic_exposed"][0]

    click.echo(
        f"Selected LYS: chain {top_lys['chain']} res {top_lys['residue_number']} "
        f"(SASA {top_lys['sasa_total']})"
    )
    click.echo(
        f"Selected {top_acid['residue_type']}: chain {top_acid['chain']} "
        f"res {top_acid['residue_number']} (SASA {top_acid['sasa_total']})"
    )

    # --- Extract sequences and residue-number → index maps -----------------
    chains1 = extract_chains(pdb1)
    chains2 = extract_chains(pdb2)

    # Assign unique Boltz chain IDs across both PDBs
    boltz_id_iter = iter(ascii_uppercase)
    pdb1_boltz: dict[str, str] = {c: next(boltz_id_iter) for c in chains1}
    pdb2_boltz: dict[str, str] = {c: next(boltz_id_iter) for c in chains2}

    # --- Sequences section -------------------------------------------------
    sequences: list[dict] = []
    for cid, (seq, _) in chains1.items():
        sequences.append({"protein": {"id": pdb1_boltz[cid], "sequence": seq}})
    for cid, (seq, _) in chains2.items():
        sequences.append({"protein": {"id": pdb2_boltz[cid], "sequence": seq}})

    # --- Bond constraint ---------------------------------------------------
    lys_chain = top_lys["chain"]
    lys_resnum = top_lys["residue_number"]
    _, lys_map = chains1[lys_chain]
    lys_idx = lys_map.get(str(lys_resnum)) or lys_map.get(int(lys_resnum))

    if lys_idx is None:
        potential_keys = [k for k in lys_map.keys() if k.endswith(str(lys_resnum))]
        if potential_keys:
            lys_idx = lys_map[potential_keys[0]]

    acid_chain = top_acid["chain"]
    acid_resnum = str(top_acid["residue_number"])
    _, acid_map = chains2[acid_chain]
    acid_idx = acid_map[acid_resnum]
    acid_atom = CARBOXYL_ATOM[top_acid["residue_type"]]

    constraints = [
        {
            "bond": {
                "atom1": [pdb1_boltz[lys_chain], lys_idx, LYS_ATOM],
                "atom2": [pdb2_boltz[acid_chain], acid_idx, acid_atom],
            }
        }
    ]

    # --- Write output ------------------------------------------------------
    boltz_input = {
        "version": 1,
        "sequences": sequences,
        "constraints": constraints,
    }

    with open(output, "w") as fh:
        yaml.dump(boltz_input, fh, default_flow_style=False, sort_keys=False)

    click.echo(f"Boltz-2 input written to {output}")


if __name__ == "__main__":
    main()
