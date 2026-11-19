"""Write data/raw/MANIFEST.json: SHA-256 + byte size of every data/catalogues/wildcards file."""
import hashlib
import json
from pathlib import Path

raw = Path("data/raw/tb-bdq-cat-main")
files = {}
for sub in ("data", "catalogues", "wildcards"):
    for p in sorted((raw / sub).rglob("*")):
        if p.is_file() and p.name != ".DS_Store":
            files[str(p.relative_to(raw)).replace("\\", "/")] = {
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                "bytes": p.stat().st_size,
            }
manifest = {
    "files": files,
    "manifest_sha256": hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
}
Path("data/raw/MANIFEST.json").write_text(json.dumps(manifest, indent=2))
print("manifest files:", len(files), "| manifest_sha256:", manifest["manifest_sha256"][:16])