from pathlib import Path
from typing import Optional

import modal

here = Path(__file__).parent  # the directory of this file

MINUTES = 60  # seconds

app = modal.App(name="example-boltz-predict")


@app.local_entrypoint()
def main(
    force_download: bool = False, input_yaml_path: Optional[str] = None, args: str = ""
):
    print("🧬 loading model remotely")
    download_model.remote(force_download)
    print(input_yaml_path)
    if input_yaml_path is None:
        input_yaml_path = here / "data" / "boltz_affinity.yaml"
    else:
        input_yaml_path = Path(input_yaml_path)
    input_yaml = input_yaml_path.read_text()

    print(f"🧬 running boltz with input from {input_yaml_path}")
    output = boltz_inference.remote(input_yaml)

    output_path = Path("/tmp") / "boltz" / "boltz_result.tar.gz" # CHANGE THIS I DON=T LIKE IT
    output_path.parent.mkdir(exist_ok=True, parents=True)
    print(f"🧬 writing output to {output_path}")
    output_path.write_bytes(output)

image = modal.Image.debian_slim(python_version="3.12").uv_pip_install("boltz==2.1.1")

boltz_model_volume = modal.Volume.from_name("boltz-models", create_if_missing=True)
models_dir = Path("/models/boltz")

# FIX: Define the uploaded templates volume 
boltz_template_volume = modal.Volume.from_name("boltz-templates")
templates_dir = Path("/templates")

@app.function(
    image=image,
    volumes={models_dir: boltz_model_volume,
    templates_dir: boltz_template_volume,},
    timeout=10 * MINUTES,
    gpu="H100",
)
def boltz_inference(boltz_input_yaml: str, args="") -> bytes:
    import shlex
    import subprocess

    input_path = Path("input.yaml")
    input_path.write_text(boltz_input_yaml)

    args = shlex.split(args)

    print(f"🧬 predicting structure using boltz model from {models_dir}")
    subprocess.run(
        ["boltz", "predict", input_path, "--use_msa_server", "--cache", str(models_dir)]
        + args,
        check=True,
    )

    print("🧬 packaging up outputs")
    output_bytes = package_outputs(f"boltz_results_{input_path.with_suffix('').name}")

    return output_bytes


def package_outputs(output_dir: str) -> bytes:
    import io
    import tarfile

    tar_buffer = io.BytesIO()

    with tarfile.open(fileobj=tar_buffer, mode="w:gz") as tar:
        tar.add(output_dir, arcname=output_dir)

    return tar_buffer.getvalue()

download_image = (
    modal.Image.debian_slim()
    .uv_pip_install("huggingface-hub==0.36.0")
    .env({"HF_XET_HIGH_PERFORMANCE": "1"})
)


@app.function(
    volumes={models_dir: boltz_model_volume},
    timeout=20 * MINUTES,
    image=download_image,
)
def download_model(
    force_download: bool = False,
    revision: str = "6fdef46d763fee7fbb83ca5501ccceff43b85607",
):
    from huggingface_hub import snapshot_download

    snapshot_download(
        repo_id="boltz-community/boltz-2",
        revision=revision,
        local_dir=models_dir,
        force_download=force_download,
    )
    boltz_model_volume.commit()

    print(f"🧬 model downloaded to {models_dir}")


