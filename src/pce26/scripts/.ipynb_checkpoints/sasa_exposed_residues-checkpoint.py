#!/usr/bin/env python3
"""Compute exposed acidic (GLU/ASP) and basic (LYS) residues from two PDB files using FreeSASA."""

from pathlib import Path

import click
import freesasa
import MDAnalysis as mda
import yaml


ACIDIC_RESIDUES = {"GLU", "ASP"}
BASIC_RESIDUES = {"LYS"}
TARGET_RESIDUES = ACIDIC_RESIDUES | BASIC_RESIDUES

DEFAULT_THRESHOLD = 0.2  # 20% relative side-chain SASA


def compute_exposed_residues(
    pdb_path: str, threshold: float
) -> dict[str, list[dict]]:
    """Return exposed GLU/ASP and LYS residues sorted by SASA (descending)."""
    u = mda.Universe(pdb_path)
    structure = freesasa.Structure()
    for atom in u.atoms:
        structure.addAtom(
            atom.name,
            atom.resname,
            str(atom.resid),
            atom.chainID,
            float(atom.position[0]),
            float(atom.position[1]),
            float(atom.position[2]),
        )
    result = freesasa.calc(structure)
    residue_areas = result.residueAreas()

    acidic: list[dict] = []
    basic: list[dict] = []

    for chain_label, residues in residue_areas.items():
        for res_number, area in residues.items():
            if area.residueType not in TARGET_RESIDUES:
                continue
            if not area.hasRelativeAreas:
                continue
            if area.relativeSideChain < threshold:
                continue

            entry = {
                "chain": chain_label,
                "residue_number": res_number,
                "residue_type": area.residueType,
                "sasa_total": round(area.total, 2),
                "sasa_side_chain": round(area.sideChain, 2),
                "relative_side_chain": round(area.relativeSideChain, 3),
            }

            if area.residueType in ACIDIC_RESIDUES:
                acidic.append(entry)
            else:
                basic.append(entry)

    # Rank by total SASA descending
    acidic.sort(key=lambda r: r["sasa_total"], reverse=True)
    basic.sort(key=lambda r: r["sasa_total"], reverse=True)

    return {"acidic_exposed": acidic, "basic_exposed": basic}


@click.command()
@click.argument("input1", type=click.Path(exists=True, dir_okay=False))
@click.argument("input2", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o",
    "--output",
    default="exposed_residues.yaml",
    show_default=True,
    help="Output YAML file path.",
)
@click.option(
    "-t",
    "--threshold",
    default=DEFAULT_THRESHOLD,
    show_default=True,
    type=float,
    help="Minimum relative side-chain SASA to consider a residue exposed (0-1).",
)
def main(input1: str, input2: str, output: str, threshold: float) -> None:
    """Identify exposed GLU/ASP and LYS residues in two PDB files via FreeSASA."""
    results: dict = {"structures": []}

    for pdb_path in (input1, input2):
        name = Path(pdb_path).name
        click.echo(f"Processing {name} ...")
        exposed = compute_exposed_residues(pdb_path, threshold)
        results["structures"].append({"pdb": name, **exposed})

    with open(output, "w") as fh:
        yaml.dump(results, fh, default_flow_style=False, sort_keys=False)

    click.echo(f"Results written to {output}")
   
    #To specify the top exposed residues for each structures
    
    str1 = results["structures"][0]
    str2 = results["structures"][1]
    
    
    if str1["basic_exposed"]:
        top_lys = str1["basic_exposed"][0]
        click.echo(
            f"{str1['pdb']} - most exposed lysine:\n"
            f"- Chain: {top_lys['chain']}\n"
            f"- Residue index: {top_lys['residue_number']}\n"
            f"- Total SASA: {top_lys['sasa_total']}"
         )
    else:
        click.echo(f"{str1['pdb']} - no exposed lysine found.")

    if str2["acidic_exposed"]:
        top_acid = str2["acidic_exposed"][0]
        click.echo(
            f"{str2['pdb']} - most exposed acidic residue:\n"
            f"- Chain: {top_acid['chain']} \n"
            f"- Residue index: {top_acid['residue_number']} \n"
            f"- Total SASA: {top_acid['sasa_total']}"
        )
    else:
        click.echo(f"{str2['pdb']} - no exposed GLU/ASP found.")

if __name__ == "__main__":
    main()
