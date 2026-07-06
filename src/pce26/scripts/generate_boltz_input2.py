import numpy as np
import pandas as pd
from Bio import PDB
import yaml
import matplotlib.pyplot as plt
from string import ascii_uppercase
import MDAnalysis as mda
import argparse
import sys
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

def kmeans(yaml_path, assembled_pdb_path):
    # Load SASA result and create a sasa_map
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)

        # Extract the list of lysines
        lys = data['structures'][0]['basic_exposed']

        # Create a dictionary with chain_id, residue_number, and sasa_total
        sasa_map = {}
        for item in lys:
            chain = item['chain']
            residue_num = item['residue_number']
            sasa_map[(chain, residue_num)] = item['sasa_total'] 
    
            print(f"Loaded {len(sasa_map)} residues with SASA scores.")

    # Parse PDB 
    parser = PDB.PDBParser()

    # Load the structure
    structure = parser.get_structure("The assembled virus", assembled_pdb_path)

    print(f"Structure {structure.id} loaded successfully!")
      
    final_features = [] # features for KMeans

    for model in structure:
        for chain in model:
            c_id = chain.id.strip().upper()

            for residue in chain:
                if residue.get_resname() == "LYS":
                    res_num = residue.get_id()[1]

                    # sasa residue_number is str, and res_num from pdb file is int. 
                    if (c_id, str(res_num)) in sasa_map:
                        sasa_val = sasa_map[(c_id, str(res_num))]
                    
                        # Calculate distances
                        coords = residue['CA'].get_coord()
                        dist = (coords[0]**2 + coords[1]**2 + coords[2]**2)**0.5
                    
                        final_features.append([c_id, res_num, dist, sasa_val])
               
                        print(f"Lysine at {res_num} is {dist:.2f} Angstroms from center")

    # Convert to a final array for KMeans
    X_combined = np.array(final_features)
    X_combined

    # KMeans input
    X_input = X_combined[:, 2:].astype(float) 

    # Scale the distance and sasa scores
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_input)

    # Run KMeans
    kmeans = KMeans(n_clusters=2, random_state=1234, n_init="auto")
    kmeans.fit(X_scaled)

    # Labels the original data
    labels = kmeans.labels_

    #Results & Visualization
    df = pd.DataFrame(X_combined, columns=['Chain', 'ResNum', 'Dist', 'SASA'])
    df['Cluster'] = labels
    df[['Dist', 'SASA']] = df[['Dist', 'SASA']].astype(float)

    # Print Cluster Summaries
    for i in range(2):
        mask = (labels == i)
        print(f"Cluster {i}: Avg Dist = {X_input[mask,0].mean():.2f} Å, Avg SASA = {X_input[mask,1].mean():.2f}")

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.scatter(X_input[:, 0], X_input[:, 1], c=labels, cmap='viridis', alpha=0.6)
    plt.colorbar(label='Cluster ID')
    plt.xlabel('Distance from Center (Å)')
    plt.ylabel('SASA Score')
    plt.title(f'Lysine Clustering: {assembled_pdb_path.split("/")[-1]}')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()

    # Find the lysine with biggest sasa score in the cluster with higher average distance

    # Calculate the average distance for both clusters
    dist_0 = X_input[labels == 0, 0].mean()
    dist_1 = X_input[labels == 1, 0].mean()

    # Identify the outward cluster label 
    if dist_1 > dist_0:
        outward_label = 1
    else:
        outward_label = 0

    print(f"Dynamic check: Cluster {outward_label} is the Outward shell.")

    # Create X_combined for cluster 1
    X_combined_c1 = X_combined[labels == outward_label]

    # Cluster 1
    cluster_1 = X_input[labels==outward_label,]
    cluster_1

    # Find the one with the highest sasa score
    max_idx = np.argmax(cluster_1[:, 1])

    top_lys = X_combined_c1[max_idx]
    top_lys_chain = top_lys[0]
    top_lys_res = top_lys[1]
    top_lys_rest = "NZ"

    print(f"Top Lysine for Boltz: Chain {top_lys_chain}, Residue {top_lys_res}")
    
    # Choosing acid which shows the highest sasa score
    # Acidic data 
    acidic =  data['structures'][1]['acidic_exposed']

    # Choose the one with the highest sasa score
    top_acid = acidic[0]
    top_acid
    top_acid_chain = top_acid['chain']
    top_acid_res = top_acid['residue_number']

    top_acid_rest=" "
    if top_acid['residue_type'] == 'GLU':
        top_acid_rest = "CD"
    elif top_acid['residue_type'] == 'ASP':
        top_acid_rest = "CG"

    return {
        "lys_chain": top_lys_chain,
        "lys_res": top_lys_res,
        "lys_atom": "NZ",
        "acid_chain": top_acid['chain'],
        "acid_res": str(top_acid['residue_number']),
        "acid_atom": top_acid_rest
    }

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
    parser = argparse.ArgumentParser(description='Generate Boltz-2 input from SASA Clustering')
    parser.add_argument('--yaml', required=True, help='Path to SASA YAML')
    parser.add_argument('--asmbl', required=True, help='Path to the assembled PDB')
    parser.add_argument('--pdb1', required=True, help='Path to first PDB (Lysine source)')
    parser.add_argument('--pdb2', required=True, help='Path to second PDB (Acidic source)')
    parser.add_argument('--output', default='boltz_input.yaml', help='Output file name')
    args = parser.parse_args()

    # Run Analysis
    targets = kmeans(args.yaml, args.asmbl)

    # Extract Sequences
    chains1 = extract_chains(args.pdb1)
    chains2 = extract_chains(args.pdb2)

    # Assign IDs
    boltz_id_iter = iter(ascii_uppercase)
    pdb1_mapping = {c: next(boltz_id_iter) for c in chains1}
    pdb2_mapping = {c: next(boltz_id_iter) for c in chains2}

    # uild Sequence List
    sequences = []
    for cid, (seq, _) in chains1.items():
        sequences.append({"protein": {"id": pdb1_mapping[cid], "sequence": seq}})
    for cid, (seq, _) in chains2.items():
        sequences.append({"protein": {"id": pdb2_mapping[cid], "sequence": seq}})

    # Build Constraints
    # Get 1-based index from MDAnalysis map
    try:
        lys_idx = chains1[targets['lys_chain']][1][targets['lys_res']]
        acid_idx = chains2[targets['acid_chain']][1][targets['acid_res']]
    except KeyError as e:
        print(f"Error: Residue mapping failed. Check if PDB and YAML numbering match. Missing: {e}")
        sys.exit(1)

    constraints = [{
        "bond": {
            "atom1": [pdb1_mapping[targets['lys_chain']], lys_idx, targets['lys_atom']],
            "atom2": [pdb2_mapping[targets['acid_chain']], acid_idx, targets['acid_atom']],
        }
    }]

    # Export
    boltz_input = {
        "version": 1,
        "sequences": sequences,
        "constraints": constraints,
    }

    with open(args.output, "w") as f:
        yaml.dump(boltz_input, f, default_flow_style=False, sort_keys=False)

    print(f"Successfully generated {args.output}")
    print(f"Bond: {targets['lys_chain']}:{targets['lys_res']} (Lys) -> {targets['acid_chain']}:{targets['acid_res']} (Acid)")

if __name__ == "__main__":
    main()