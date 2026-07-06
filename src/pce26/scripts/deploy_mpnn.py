from pathlib import Path
from typing import Optional

import modal

MINUTES = 60  # seconds

app = modal.App(name="protein-mpnn")

# The official Rosetta Commons ProteinMPNN image already ships with the model
# weights baked in (at /app/proteinmpnn/vanilla_model_weights, plus the
# CA-only and soluble variants). There is therefore no need for a separate
# weight-download step or a Modal Volume to persist them, unlike Boltz-2.
#
# We clear the inherited ENTRYPOINT (which is `python /app/proteinmpnn/
# protein_mpnn_run.py`) so that Modal's container entrypoint can run
# instead; otherwise Modal's CLI arguments get appended to ProteinMPNN's
# argv and the run fails with "unrecognized arguments".
image = (
    modal.Image.from_registry("rosettacommons/proteinmpnn")
    .dockerfile_commands(["ENTRYPOINT []"])
)

# Path of the ProteinMPNN source tree inside the container.
MPNN_DIR = "/app/proteinmpnn"


@app.function(
    image=image,
    timeout=30 * MINUTES,
    gpu="A10G",
)
def mpnn_inference(
    pdb_text: str,
    chains: Optional[str] = None,
    args: str = "",
) -> bytes:
    """Run ProteinMPNN on a single PDB file passed as text.

    Parameters
    ----------
    pdb_text : str
        Raw contents of the input PDB file.
    chains : Optional[str]
        Space-separated list of chain IDs to redesign, e.g. ``"A B"``. If
        ``None``, all chains in the PDB are designed.
    args : str
        Extra CLI flags forwarded verbatim to ``protein_mpnn_run.py`` (e.g.
        ``"--num_seq_per_target 8 --sampling_temp 0.1 --model_name v_48_020"``).
    """
    import shlex
    import subprocess

    work_dir = Path("/tmp/mpnn_run")
    work_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = work_dir / "input.pdb"
    pdb_path.write_text(pdb_text)

    out_folder = work_dir / "output"
    out_folder.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python",
        f"{MPNN_DIR}/protein_mpnn_run.py",
        "--pdb_path", str(pdb_path),
        "--out_folder", str(out_folder),
    ]
    if chains:
        cmd += ["--pdb_path_chains", chains]
    cmd += shlex.split(args)

    print(f"🧬 running ProteinMPNN: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    print("🧬 packaging up outputs")
    return _package_outputs(out_folder)


def _package_outputs(output_dir: Path) -> bytes:
    import io
    import tarfile

    tar_buffer = io.BytesIO()
    with tarfile.open(fileobj=tar_buffer, mode="w:gz") as tar:
        tar.add(str(output_dir), arcname=output_dir.name)
    return tar_buffer.getvalue()


@app.local_entrypoint()
def main(
    input_path: str,
    chains: Optional[str] = None,
    output_path: Optional[str] = None,
    args: str = "",
):
    """Run ProteinMPNN on Modal.

    Example
    -------
    modal run deploy_mpnn_v2.py \
        --input-path my_protein.pdb \
        --chains "A B" \
        --args "--num_seq_per_target 8 --sampling_temp 0.1"
    """
    pdb_path = Path(input_path)
    pdb_text = pdb_path.read_text()

    print(
        f"🧬 running ProteinMPNN on {pdb_path}"
        + (f" (designing chains: {chains})" if chains else " (designing all chains)")
    )
    output = mpnn_inference.remote(pdb_text, chains=chains, args=args)

    if output_path is None:
        out = Path("/tmp") / "protein-mpnn" / f"{pdb_path.stem}_mpnn.tar.gz"
    else:
        out = Path(output_path)
    out.parent.mkdir(exist_ok=True, parents=True)
    print(f"🧬 writing output to {out}")
    out.write_bytes(output)
