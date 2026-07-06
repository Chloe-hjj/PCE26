#!/usr/bin/env python3
"""Build an icosahedral assembly from a reference PDB and a designed protein.

Reads BIOMT symmetry operators from a *reference* PDB (via ``gemmi``), aligns a
*designed* protein onto each chain of the reference asymmetric unit (ASU) using
Cα RMSD superposition (via ``MDAnalysis``), then replicates the aligned ASU with
the BIOMT transforms to produce the full particle.

Typical usage for a T=3 capsid (e.g. 7MRY)::

    build-icosahedral 7MRY.pdb designed.pdb -o assembly.pdb

This aligns ``designed.pdb`` to each of the three ASU chains (A, B, C) and
applies the 60 BIOMT operators → 180-subunit particle.

Limitations
-----------
* PDB format supports at most 62 unique single-character chain IDs
  (A–Z, a–z, 0–9).  When the assembly exceeds this, chain IDs are recycled
  and the four-character **segment ID** (``segid``) distinguishes symmetry
  copies (``C00``–``C59``).
* Atom serial numbers wrap at 99 999 (PDB format limit).
"""

import string
from typing import Sequence

import click
import gemmi
import MDAnalysis as mda
import numpy as np
from MDAnalysis import Merge
from MDAnalysis.analysis import align

# Single-character chain IDs available in PDB format.
_CHAIN_IDS = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)


# ---------------------------------------------------------------------------
# BIOMT helpers
# ---------------------------------------------------------------------------

def parse_biomt(
    pdb_path: str,
) -> tuple[list[tuple[np.ndarray, np.ndarray]], list[str]]:
    """Extract BIOMT operators and target chain IDs from a PDB file.

    Returns
    -------
    operators
        List of *(R, t)* tuples — 3×3 rotation matrix and 3-vector translation.
    chains
        Chain IDs the operators should be applied to.
    """
    st = gemmi.read_structure(pdb_path)
    if not st.assemblies:
        raise click.ClickException(
            f"No biological assembly (REMARK 350) found in {pdb_path}."
        )
    assembly = st.assemblies[0]

    operators: list[tuple[np.ndarray, np.ndarray]] = []
    chains: list[str] = []
   
   generators = []
    for gen in assembly.generators:
        ops = 
        chains.extend(gen.chains)
        for op in gen.operators:
            R = np.array(op.transform.mat.tolist())
            t = np.array(op.transform.vec.tolist())
            operators.append((R, t))

    return operators, sorted(set(chains))


# ---------------------------------------------------------------------------
# Alignment helper
# ---------------------------------------------------------------------------

def align_to_chains(
    ref_universe: mda.Universe,
    design_path: str,
    chain_ids: Sequence[str],
) -> list[np.ndarray]:
    """Align the designed protein to each *chain_id* of the reference.

    For every requested chain the function:
    1. Loads a fresh copy of the designed structure.
    2. Superimposes it (Cα RMSD) onto the reference chain.
    3. Stores the resulting coordinates.

    Returns a list of coordinate arrays, one per chain.
    """
    positioned: list[np.ndarray] = []
    for cid in chain_ids:
        ref_ca = ref_universe.select_atoms(f"chainID {cid} and name CA")
        if len(ref_ca) == 0:
            raise click.ClickException(
                f"No Cα atoms found for chain {cid} in reference."
            )

        mobile = mda.Universe(design_path)
        mob_ca = mobile.select_atoms("name CA")

        # Match Cα atoms by sequential position (not residue ID) so that
        # chains with different numbering schemes or a few missing
        # residues can still be aligned.
        n_common = min(len(ref_ca), len(mob_ca))
        ref_match = ref_universe.atoms[ref_ca.indices[:n_common]]
        mob_match = mobile.atoms[mob_ca.indices[:n_common]]

        _, rmsd = align.alignto(mob_match, ref_match, select="all")

        n_dropped = max(len(ref_ca), len(mob_ca)) - n_common
        extra = f" ({n_dropped} residues excluded)" if n_dropped else ""
        click.echo(
            f"  Chain {cid}: RMSD = {rmsd:.2f} Å "
            f"({n_common} Cα pairs{extra})"
        )
        positioned.append(mobile.atoms.positions.copy())

    return positioned


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command()
@click.argument("reference_pdb", type=click.Path(exists=True, dir_okay=False))
@click.argument("designed_pdb", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "-o",
    "--output",
    default="assembly.pdb",
    show_default=True,
    help="Output PDB file path.",
)
@click.option(
    "-c",
    "--chains",
    default=None,
    help="Comma-separated reference chain IDs to align against (default: all "
    "chains listed in the BIOMT records).",
)
@click.option(
    "-n",
    "--n-copies",
    default=0,
    show_default=False,
    type=int,
    help="Limit the number of BIOMT operators applied (0 = all).",
)
def main(
    reference_pdb: str,
    designed_pdb: str,
    output: str,
    chains: str | None,
    n_copies: int,
) -> None:
    """Build an icosahedral assembly from REFERENCE_PDB symmetry and DESIGNED_PDB.

    REFERENCE_PDB provides the BIOMT symmetry operators and the alignment
    target (asymmetric unit).  DESIGNED_PDB is structurally aligned onto each
    requested chain of the reference ASU, then replicated.
    """
    # -- 1. Parse BIOMT operators from the reference -------------------------
    operators, biomt_chains = parse_biomt(reference_pdb)
    if n_copies > 0:
        operators = operators[:n_copies]
    click.echo(
        f"Reference: {len(operators)} BIOMT operators, "
        f"chains {', '.join(biomt_chains)}"
    )

    # -- 2. Determine which chains to align against -------------------------
    if chains is not None:
        align_chains = [c.strip() for c in chains.split(",")]
    else:
        align_chains = biomt_chains
    click.echo(f"Aligning designed protein to chain(s): {', '.join(align_chains)}")

    # -- 3. Align designed protein to each reference chain ------------------
    ref_u = mda.Universe(reference_pdb)
    positioned_coords = align_to_chains(ref_u, designed_pdb, align_chains)

    # -- 4. Build a template Universe for one ASU copy ----------------------
    #    Load designed protein once per chain and merge into a single ASU.
    asu_parts: list[mda.AtomGroup] = []
    for i, coords in enumerate(positioned_coords):
        part = mda.Universe(designed_pdb)
        part.atoms.positions = coords
        # Pre-assign chain IDs so the ASU chains are distinguishable
        for seg in part.segments:
            seg.atoms.chainIDs = np.array(
                [align_chains[i]] * len(seg.atoms)
            )
        asu_parts.append(part.atoms)

    asu_merged = Merge(*asu_parts)
    asu_coords = asu_merged.atoms.positions.copy()
    n_asu_atoms = len(asu_coords)
    n_asu_chains = len(align_chains)

    click.echo(
        f"Designed ASU: {n_asu_atoms} atoms, {n_asu_chains} chain(s)"
    )

    # -- 5. Replicate using BIOMT transforms --------------------------------
    click.echo(f"Applying {len(operators)} BIOMT symmetry operations ...")

    copies: list[mda.AtomGroup] = []
    
    for i, (R, t) in enumerate(operators, start=1):
        copy = asu_merged.copy()
        copy.atoms.positions = (R @ asu_coords.T).T + t

        for seg in copy.segments:
            #cid = _CHAIN_IDS[chain_counter % len(_CHAIN_IDS)]
            #seg.atoms.chainIDs = np.array([cid] * len(seg.atoms))
            #seg.segid = f"C{i:02d}"
            #chain_counter += 1
            original_chain = seg.atoms.chainIDs[0]
            seg.atoms.chainIDs = np.array([original_chain]*len(seg.atoms))
            seg.segid = f"{original_chain}{i}"

        copies.append(copy.atoms)

    merged = Merge(*copies)
    n_total = len(merged.atoms)
    total_chains = n_asu_chains * len(operators)
    click.echo(
        f"Assembly: {n_asu_atoms} atoms × {len(operators)} copies "
        f"= {n_total} atoms, {total_chains} chain(s)"
    )

    merged.atoms.write(output)
    click.echo(f"Assembly written to {output}")


if __name__ == "__main__":
    main()
