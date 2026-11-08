"""Geological-origin audit and LGM-connectivity crosswalk for fd.

Uses independent, published Weigelt 2013 island-centroid/area/GMMC data
and the Zizka 2022 supplementary *missing-geology additions only*.
Neither great-circle distance nor GMMC=0 is geological proof of oceanic origin.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import urllib.request
import zipfile

import numpy as np

RADIUS_KM = 6371.0088
WEIGELT = "https://zenodo.org/records/4971536/files/Weigelt_etal_2013_PNAS_islanddata.csv?download=1"
ZIZKA = "https://zenodo.org/records/6325640/files/input.zip?download=1"
GSHHG_COMMIT = "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c"
GSHHG_URL = ("https://raw.githubusercontent.com/zuizui0223/island/" + GSHHG_COMMIT
             + "/results/geography_20260924/corrected_geography_covariates.csv")


def download(url: str, dest: Path, max_bytes: int) -> str:
    sha = hashlib.sha256()
    size = 0
    request = urllib.request.Request(url, headers={"User-Agent": "fd-geology-audit/1.0"})
    with urllib.request.urlopen(request, timeout=120) as source, dest.open("wb") as out:
        while chunk := source.read(1 << 20):
            size += len(chunk)
            if size > max_bytes:
                raise ValueError("Source exceeds declared maximum byte count")
            sha.update(chunk)
            out.write(chunk)
    if size == 0:
        raise ValueError("Downloaded empty source")
    return sha.hexdigest()


def _number(value) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def _csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle)


def xyz(lat, lon) -> np.ndarray:
    a, b = np.radians(np.asarray(lat, dtype=float)), np.radians(np.asarray(lon, dtype=float))
    return np.column_stack((np.cos(a) * np.cos(b), np.cos(a) * np.sin(b), np.sin(a)))


def _valid_features(lat: float, lon: float, area: float) -> bool:
    return bool(np.isfinite([lat, lon, area]).all() and -90 <= lat <= 90
                and -180 <= lon <= 180 and area > 0)


def parse_gshhg(path: Path) -> list[dict]:
    rows = []
    for r in _csv(path):
        lat, lon, area = (_number(r.get(k)) for k in
                          ("island_latitude", "island_longitude", "area_km2"))
        if _valid_features(lat, lon, area):
            rows.append({"id": r["island_id"], "lat": lat, "lon": lon, "area": area,
                         "regime": r.get("analysis_regime", ""),
                         "continent_distance_km": _number(r.get("distance_to_continent_km"))})
    if not rows or len({r["id"] for r in rows}) != len(rows):
        raise ValueError("No valid unique geographic islands")
    return rows


def parse_weigelt(path: Path) -> list[dict]:
    rows = []
    for r in _csv(path):
        lat, lon, area = (_number(r.get(k)) for k in ("Lat", "Long", "Area"))
        if _valid_features(lat, lon, area) and r.get("GMMC", "") in ("0", "1"):
            rows.append({"id": str(r["ID"]), "lat": lat, "lon": lon, "area": area,
                         "connection": "lgm_connected" if r["GMMC"] == "1"
                         else "lgm_unconnected", "name": r.get("Island", ""),
                         "archipelago": r.get("Archip", "")})
    return rows


def parse_dbf(buf: bytes) -> list[dict]:
    if len(buf) < 32:
        raise ValueError("Invalid DBF")
    n = int.from_bytes(buf[4:8], "little")
    header = int.from_bytes(buf[8:10], "little")
    record = int.from_bytes(buf[10:12], "little")
    fields = []
    offset = 1
    for position in range(32, header - 1, 32):
        desc = buf[position:position + 32]
        if not desc or desc[0] == 13:
            break
        name = desc[:11].split(bytes([0]), 1)[0].decode("ascii", errors="replace")
        fields.append((name, offset, desc[16]))
        offset += desc[16]
    if offset > record:
        raise ValueError("DBF field offsets exceed record size")
    result = []
    for i in range(n):
        row_bytes = buf[header + i * record:header + (i + 1) * record]
        if len(row_bytes) != record:
            raise ValueError("Truncated DBF rows")
        if row_bytes[0:1] == b"*":
            continue
        row = {name: row_bytes[start:start + width].decode("latin1").strip()
               for name, start, width in fields}
        result.append(row)
    return result


def parse_verified_geology(archive: Path) -> list[dict]:
    with zipfile.ZipFile(archive) as z:
        with z.open("input/entities_with_missing_geology_added.csv") as handle:
            manual = list(csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8-sig")))
        entities = {r["entt_ID"]: r for r in parse_dbf(z.read("input/geoentities_simple.dbf"))}
    rows = []
    for r in manual:
        entity = entities.get(r["entity_ID"])
        if not entity:
            continue
        label = r.get("geology_simple_add", "").strip().lower()
        # Preserve the geological distinction, with mixed and unknown not
        # automatically forced into a binary oceanic/continental class.
        if label in ("oceanic", "volcanic", "atoll", "raised ocean floor"):
            kind = "oceanic"
        elif label in ("continental", "continental shelf", "shelf", "continental fragment"):
            kind = "continental"
        elif label == "mixed":
            kind = "mixed"
        else:
            kind = "unresolved"
        lat, lon, area = (_number(entity.get(k)) for k in ("point_y", "point_x", "area"))
        if _valid_features(lat, lon, area):
            rows.append({"id": r["entity_ID"], "lat": lat, "lon": lon, "area": area,
                         "geological_origin": kind,
                         "geological_source_value": label,
                         "geological_source_name": r.get("geo_entity", ""),
                         "geological_reference": r.get("reference", "")})
    return rows


def crosswalk(target: list[dict], source: list[dict], *, k: int = 24) -> tuple[dict, dict]:
    """Conservative many-to-one ambiguity-screened spherical centroid crosswalk.

    A candidate must match area within a factor 1.8 and centroid distance
    under a predeclared size-adaptive bound. Reused source IDs rejected.
    """
    if not source:
        return {}, {"matched": 0, "reason": "empty source"}
    features = xyz([r["lat"] for r in source], [r["lon"] for r in source])
    result = {}
    uncertain = {"no_candidate": 0, "ambiguous": 0, "source_reused": 0}
    for begin in range(0, len(target), 80):
        subset = target[begin: begin + 80]
        query = xyz([r["lat"] for r in subset], [r["lon"] for r in subset])
        similarities = query @ features.T
        topn = min(k, len(source))
        indices = np.argpartition(-similarities, topn - 1, axis=1)[:, :topn]
        for local, row in enumerate(subset):
            candidates = []
            limit = max(9.0, min(70.0, 0.32 * math.sqrt(row["area"])))
            for index in indices[local]:
                other = source[int(index)]
                d = RADIUS_KM * math.acos(float(np.clip(similarities[local, index], -1.0, 1.0)))
                ratio = max(row["area"], other["area"]) / min(row["area"], other["area"])
                if d <= limit and ratio <= 1.8:
                    score = d / limit + abs(math.log(row["area"] / other["area"])) / math.log(1.8)
                    candidates.append((score, d, ratio, other))
            if not candidates:
                uncertain["no_candidate"] += 1
                continue
            candidates.sort(key=lambda x: (x[0], x[1], x[3]["id"]))
            best = candidates[0]
            # A second plausible island is problematic even if the best has
            # a slightly better score; do not force an arbitrary name join.
            if len(candidates) > 1 and candidates[1][0] - best[0] < 0.20:
                uncertain["ambiguous"] += 1
                continue
            result[row["id"]] = {"source": best[3], "km": best[1],
                                 "area_ratio": best[2], "match_score": best[0],
                                 "candidate_count": len(candidates)}
    reverse = {}
    for island_id, info in result.items():
        reverse.setdefault(info["source"]["id"], []).append(island_id)
    for ids in reverse.values():
        if len(ids) > 1:
            for island_id in ids:
                del result[island_id]
                uncertain["source_reused"] += 1
    audit = {"target_total": len(target), "source_total": len(source),
             "matched": len(result), **uncertain,
             "max_area_ratio": 1.8, "min_score_gap_second_best": 0.20,
             "max_centroid_km": "max(9,min(70,0.32*sqrt(area_km2)))",
             "note": "conservative centroid+area join, NOT verified polygon identity"}
    return result, audit


def run_crosswalk(destination: Path) -> dict:
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as work:
        root = Path(work)
        files = {"gshhg": root / "gshhg.csv", "weigelt": root / "weigelt.csv",
                 "zizka": root / "zizka.zip"}
        hashes = {"gshhg": download(GSHHG_URL, files["gshhg"], 80_000_000),
                  "weigelt": download(WEIGELT, files["weigelt"], 15_000_000),
                  "zizka": download(ZIZKA, files["zizka"], 100_000_000)}
        geo = parse_gshhg(files["gshhg"])
        conn = parse_weigelt(files["weigelt"])
        verified = parse_verified_geology(files["zizka"])
    conn_match, conn_audit = crosswalk(geo, conn)
    geo_match, geo_audit = crosswalk(geo, verified)
    source = {r["id"]: r for r in geo}
    rows = []
    for island in geo:
        matched = conn_match.get(island["id"])
        strict = geo_match.get(island["id"])
        direct = strict["source"] if strict else {}
        connected = matched["source"] if matched else {}
        label = direct.get("geological_origin", "unresolved")
        if label not in ("continental", "oceanic", "mixed"):
            label = "unresolved"
        rows.append({
            "island_id": island["id"],
            "analysis_regime": island["regime"],
            "area_km2": island["area"],
            "distance_to_continent_km": island["continent_distance_km"],
            "island_latitude": island["lat"], "island_longitude": island["lon"],
            "geological_origin": label,
            "geological_class_source": "Zizka_2022_missing_geology_additions" if strict else "",
            "geological_record_id": direct.get("id", ""),
            "geological_reference": direct.get("geological_reference", ""),
            "geological_match_km": strict["km"] if strict else "",
            "geological_area_ratio": strict["area_ratio"] if strict else "",
            "LGM_connectivity": connected.get("connection", "unresolved"),
            "weigelt_record_id": connected.get("id", ""),
            "weigelt_island": connected.get("name", ""),
            "weigelt_archipelago": connected.get("archipelago", ""),
            "weigelt_match_km": matched["km"] if matched else "",
            "weigelt_area_ratio": matched["area_ratio"] if matched else "",
        })
    columns = list(rows[0])
    with (destination / "island_origin_crosswalk.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    from collections import Counter
    geology_counts = dict(Counter(r["geological_origin"] for r in rows))
    connectivity_counts = dict(Counter(r["LGM_connectivity"] for r in rows))
    contradictions = [r["island_id"] for r in rows if r["geological_origin"] == "oceanic"
                      and r["LGM_connectivity"] == "lgm_connected"]
    audit = {
        "status": "geological_origin_partial_and_lgm_proxy_separate",
        "target_gshhg": len(geo), "weigelt_source_records": len(conn),
        "zizka_verified_additional_records": len(verified),
        "geological_counts": geology_counts, "LGM_counts": connectivity_counts,
        "geological_match_audit": geo_audit, "LGM_match_audit": conn_audit,
        "origin_LGM_inconsistency_count": len(contradictions),
        "origin_LGM_inconsistency_ids": contradictions[:20],
        "source_sha256": hashes,
        "interpretation_limit": "LGM-unconnected does not imply geologically oceanic",
        "source_links": {"gshhg": GSHHG_URL, "weigelt": WEIGELT, "zizka": ZIZKA},
    }
    (destination / "crosswalk_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"status": audit["status"], "geological_counts": geology_counts,
                      "LGM_counts": connectivity_counts,
                      "geology_matched": len(geo_match),
                      "LGM_matched": len(conn_match),
                      "conflicting": len(contradictions)}, ensure_ascii=False))
    return audit
