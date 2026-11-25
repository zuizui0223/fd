"""Fail closed on unreviewed active files and known own-DB input markers."""
import json,re,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/"config/external_db_ownership.json").read_text())
tracked=set(subprocess.check_output(["git","-C",str(root),"ls-files","-z"]).decode().split("\0"))-{""}
allowed=set(m["allowed_paths"])
if tracked!=allowed:
    raise SystemExit("ownership manifest mismatch: "+str({"unexpected":sorted(tracked-allowed),"missing":sorted(allowed-tracked)}))
markers=re.compile(r"chapter1_status_flora|island_species_occurrences|direct_only_trait_ledger|chapter1-progressive-analysis|floral_space_20261106|fetch_island_aggregates|get_public_flora|fd-own-db|6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",re.I)
for rel in sorted(tracked):
    if rel=="scripts/check_external_source_boundary.py":continue
    if rel.startswith(("scripts/","src/",".github/workflows/")):
        content=(root/rel).read_text()
        if markers.search(content):
            raise SystemExit("own-DB marker in active external path: "+rel)
print("External-only file allowlist and own-DB marker screen passed",len(tracked))
