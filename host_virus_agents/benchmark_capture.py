"""Append-only diagnostics for correctness runs, not a production batch store."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
from importlib import metadata


def digest(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def begin(input_file, output_file):
    sidecar = Path(str(output_file) + ".diagnostics.jsonl")
    if Path(output_file).exists() or sidecar.exists():
        raise FileExistsError("Use a new result path; existing benchmark artifacts are preserved.")
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    manifest = dict(kind="manifest", timestamp=datetime.now(timezone.utc).isoformat(),
                    input_sha256=digest(input_file), git_commit=revision,
                    python=platform.python_version(), platform=platform.platform(),
                    model_path=os.environ.get("MODEL_PATH", "./models/Qwen3-8B"),
                    sources={str(p): digest(p) for p in [*Path('.').glob('*.py'), *Path('bioresearch_env').glob('*.py')]},
                    caches_before={str(p): digest(p) for p in Path('.').glob('*cache*.json')})
    manifest["packages"] = {}
    for package in ("torch", "transformers", "requests", "tokenizers"):
        try:
            manifest["packages"][package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            manifest["packages"][package] = None
    model_dir = Path(manifest["model_path"])
    manifest["model_metadata_sha256"] = {str(p): digest(p) for p in model_dir.glob("*.json")}
    weights = [*model_dir.glob("*.safetensors"), *model_dir.glob("*.bin")]
    manifest["model_weights"] = [{"path": str(p), "size_bytes": p.stat().st_size, "sha256": digest(p)} for p in weights]
    manifest["weight_content_hashes_available"] = bool(weights)
    manifest["cache_configuration"] = {k: os.environ.get(k) for k in ("HOST_ALIAS_CACHE", "VIRUS_TAXONOMY_CACHE", "BIOLOGICAL_CONTEXT_CACHE")}
    with sidecar.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(manifest, ensure_ascii=False) + "\n")
    return sidecar


def append(sidecar, row, diagnostics):
    record = dict(kind="pair", row=row, diagnostics=diagnostics,
                  timestamp=datetime.now(timezone.utc).isoformat())
    with Path(sidecar).open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
