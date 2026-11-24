#!/usr/bin/env python3
"""Public source preflight and reproducible descriptive bird pilot, no ecological fit."""
import argparse
import csv
import hashlib
import io
import json
import pathlib
import urllib.parse
import urllib.request

DAR = "8b381ff26e3d6f4da17730dda0c5ab7dbad12eed"
BASE = "https://raw.githubusercontent.com/txm676/DARs/" + DAR + "/Data/Island_datasets/True_island_datasets/"
PILOT = {
    "Hawaii": ("Baiser et al (2017) Hawaii Birds_current_noAliens.csv", ["Vestiaria_coccinea", "Himatione_sanguinea", "Palmeria_dolei"], 6, 6, "2db621e4436856317dec96bb403bca59ac8d172b"),
    "Lesser Antilles": ("Baiser et al (2017) Lesser Antilles Birds_current_noAliens.csv", ["Eulampis_holosericeus", "Orthorhyncus_cristatus", "Eulampis_jugularis", "Cyanophaia_bicolor", "Glaucis_hirsutus"], 12, 12, "21e4db9525512375a14979479587aa3139aaf850"),
    "Marianas": ("Baiser et al (2017) Marianas Birds_current_noAliens.csv", ["Myzomela_rubratra"], 15, 12, "640e4b0f9d74197ca54a19dd2a08225f6bd18f01"),
    "Cook Islands": ("Baiser et al (2017) Cook Islands Birds_current_noAliens.csv", ["Vini_kuhlii"], 15, 1, "c14c748f59eca667d0cdbf11f7c00550781cacf6"),
    "Society Islands": ("Baiser et al (2017) Society Islands Birds_current_noAliens.csv", ["Vini_peruviana"], 11, 2, "6dda0b02a94c75b7c24b05d269a7353fb2e684a6"),
    "Wakatobi": ("OConnell et al (2020) Wakatobi.csv", ["Trichoglossus_ornatus", "Nectarinia_jugularis"], 9, 9, "a3110d9fd8a73f2dc803c655594a0a82d41058d"),
}
# Correct pinned Git blob for Wakatobi (source tree).
PILOT["Wakatobi"] = (*PILOT["Wakatobi"][:4], "a3110d9fd8a73f2dc803c655594a0a82d41058d")

def get(url, maxbytes=20_000_000):
    req = urllib.request.Request(url, headers={"User-Agent": "fd-public-source-preflight/1.0", "Accept": "*/*"})
    with urllib.request.urlopen(req, timeout=45) as r:
        data = r.read(maxbytes + 1)
    if len(data) > maxbytes:
        raise ValueError("metadata exceeds download cap")
    return data

def blobsha(b):
    return hashlib.sha1(("blob " + str(len(b)) + "\0").encode() + b).hexdigest()

def parse_matrix(raw, selected):
    rows = list(csv.reader(io.StringIO(raw.decode("utf-8-sig"))))
    if not rows or len(rows[0]) < 2 or rows[0][0].lower().strip() != "species":
        raise ValueError("unexpected matrix header")
    islands = rows[0][1:]
    if any(not x.strip() for x in islands) or len(set(islands)) != len(islands):
        raise ValueError("nonunique island header")
    species = {}
    footer = None
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        label = row[0].strip()
        if len(row) != len(islands) + 1:
            raise ValueError("unequal row size")
        if label == "Area (ha)":
            continue
        if label == "sp.r":
            if footer is not None:
                raise ValueError("duplicate footer")
            footer = [int(x) for x in row[1:]]
            continue
        if label in species or any(v.strip() not in {"0","1"} for v in row[1:]):
            raise ValueError("duplicate species or nonbinary value")
        species[label] = [int(x) for x in row[1:]]
    if footer is None or [sum(v[j] for v in species.values()) for j in range(len(islands))] != footer:
        raise ValueError("species totals disagree with source")
    if not set(selected).issubset(species):
        raise ValueError("selected taxon missing")
    counts = [sum(species[s][j] for s in selected) for j in range(len(islands))]
    return {"island_columns": len(islands), "positive_islands": sum(x > 0 for x in counts),
            "selected_taxon_count": len(selected), "selected_counts": counts,
            "selected_species_occupancy": {s: sum(species[s]) for s in selected},
            "footer_reconciled": True}

def run_pilot():
    out = {}
    for group, (path, selected, n, pos, expected_sha) in PILOT.items():
        raw = get(BASE + urllib.parse.quote(path))
        if blobsha(raw) != expected_sha:
            raise ValueError("Git blob mismatch: " + group)
        result = parse_matrix(raw, selected)
        if (result["island_columns"], result["positive_islands"]) != (n, pos):
            raise ValueError("pilot counts drift: " + group)
        result["source_file"] = path
        result["git_blob"] = expected_sha
        result["sha256"] = hashlib.sha256(raw).hexdigest()
        out[group] = result
    return out

def metadata(url):
    try:
        raw = get(url)
        return {"status": "retrieved", "url": url, "sha256": hashlib.sha256(raw).hexdigest(), "json": json.loads(raw)}
    except Exception as e:
        return {"status": "unavailable", "url": url, "error": str(e)[:400]}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.self_test:
        fake = b"species,A,B\nGenus_one,1,0\nGenus_two,1,1\nArea (ha),100,200\nsp.r,2,1\n"
        assert parse_matrix(fake, ["Genus_one"])["positive_islands"] == 1
        for bad in [fake.replace(b"sp.r,2,1", b"sp.r,1,1"), fake.replace(b"1,0", b"1,")]:
            try:
                parse_matrix(bad, ["Genus_one"])
            except ValueError:
                pass
            else:
                raise AssertionError("invalid source accepted")
        return
    # Metadata only. Do not download the multi-GB butterfly archive in this job.
    manifest = {
        "bee": metadata("https://zenodo.org/api/records/23222636"),
        "butterfly": metadata("https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.bzkh189h9"),
        "diet": metadata("https://api.figshare.com/v2/collections/3306933/articles"),
        "bird_dars_commit": DAR,
        "scope": "source intake only; no cross-guild or floral association calculated"
    }
    (args.out / "source_metadata.json").write_text(json.dumps(manifest, indent=2))
    pilot = run_pilot()
    (args.out / "bird_pilot.json").write_text(json.dumps(pilot, indent=2))
    print("Actual public pilot:", sum(r["positive_islands"] for r in pilot.values()), "/", sum(r["island_columns"] for r in pilot.values()))
    print("Metadata:", {k:v["status"] for k,v in manifest.items() if isinstance(v,dict)})
    print("No ecological H1-H4 effects estimated.")

if __name__ == "__main__":
    main()
