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
from MDAnalysis.analysis.align import rotation_matrix as calc_rotation_matrix

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
    for gen in assembly.generators:
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


def align_multichain_asu(
    ref_universe: mda.Universe,
    design_path: str,
    matching_chains: Sequence[str],
) -> mda.Universe:
    """Align a multi-chain designed ASU onto the reference ASU.

    Performs a single global superposition using Cα atoms collected from all
    *matching_chains* (design chain X → reference chain X), then applies the
    resulting rotation and translation to **every** atom in the design —
    including chains that have no counterpart in the reference.

    Returns the aligned design Universe.
    """
    design = mda.Universe(design_path)

    # Collect matched Cα coordinates per chain
    mob_ca_list: list[np.ndarray] = []
    ref_ca_list: list[np.ndarray] = []
    chain_info: list[tuple[str, int, int]] = []  # (chain, n_common, n_dropped)

    for cid in matching_chains:
        ref_ca = ref_universe.select_atoms(f"chainID {cid} and name CA")
        mob_ca = design.select_atoms(f"chainID {cid} and name CA")
        if len(ref_ca) == 0:
            raise click.ClickException(
                f"No Cα atoms found for chain {cid} in reference."
            )
        if len(mob_ca) == 0:
            raise click.ClickException(
                f"No Cα atoms found for chain {cid} in design."
            )
        n_common = min(len(ref_ca), len(mob_ca))
        mob_ca_list.append(mob_ca.positions[:n_common])
        ref_ca_list.append(ref_ca.positions[:n_common])
        n_dropped = max(len(ref_ca), len(mob_ca)) - n_common
        chain_info.append((cid, n_common, n_dropped))

    mob_all = np.vstack(mob_ca_list)
    ref_all = np.vstack(ref_ca_list)

    # Centre and compute optimal rotation
    mob_com = mob_all.mean(axis=0)
    ref_com = ref_all.mean(axis=0)
    R, global_rmsd = calc_rotation_matrix(
        mob_all - mob_com, ref_all - ref_com
    )

    # Apply transform to ALL design atoms (same method as align.alignto)
    design.atoms.translate(-mob_com)
    design.atoms.rotate(R)
    design.atoms.translate(ref_com)

    # Per-chain RMSD report (computed after alignment)
    for cid, n_common, n_dropped in chain_info:
        ref_ca = ref_universe.select_atoms(f"chainID {cid} and name CA")
        mob_ca = design.select_atoms(f"chainID {cid} and name CA")
        n = min(len(ref_ca), len(mob_ca))
        d = mob_ca.positions[:n] - ref_ca.positions[:n]
        chain_rmsd = float(np.sqrt(np.mean(np.sum(d**2, axis=1))))
        extra = f" ({n_dropped} residues excluded)" if n_dropped else ""
        click.echo(
            f"  Chain {cid}: RMSD = {chain_rmsd:.2f} Å "
            f"({n_common} Cα pairs{extra})"
        )

    total_pairs = sum(n for _, n, _ in chain_info)
    click.echo(
        f"  Global RMSD = {global_rmsd:.2f} Å ({total_pairs} total Cα pairs)"
    )

    return design


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

    # -- 3. Detect multi-chain design and build ASU -------------------------
    ref_u = mda.Universe(reference_pdb)
    design_chain_ids = sorted(set(mda.Universe(designed_pdb).atoms.chainIDs))
    matching = [c for c in align_chains if c in design_chain_ids]
    multichain = len(design_chain_ids) > 1 and len(matching) == len(align_chains)

    if multichain:
        extra = [c for c in design_chain_ids if c not in matching]
        click.echo(
            f"Multi-chain design detected "
            f"(chains {', '.join(design_chain_ids)})"
        )
        click.echo(f"Aligning matched chains: {', '.join(matching)}")
        if extra:
            click.echo(f"Extra chains included in ASU: {', '.join(extra)}")

        asu_merged = align_multichain_asu(ref_u, designed_pdb, matching)
        asu_coords = asu_merged.atoms.positions.copy()
        n_asu_atoms = len(asu_coords)
        n_asu_chains = len(design_chain_ids)
    else:
        click.echo(
            f"Aligning designed protein to chain(s): "
            f"{', '.join(align_chains)}"
        )
        positioned_coords = align_to_chains(ref_u, designed_pdb, align_chains)

        # Build a template Universe for one ASU copy.
        # Load designed protein once per chain and merge into a single ASU.
        asu_parts: list[mda.AtomGroup] = []
        for i, coords in enumerate(positioned_coords):
            part = mda.Universe(designed_pdb)
            part.atoms.positions = coords
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
    chain_counter = 0

    for i, (R, t) in enumerate(operators):
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
