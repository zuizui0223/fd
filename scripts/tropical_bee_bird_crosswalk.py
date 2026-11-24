#!/usr/bin/env python3
"""Independent native-bee / bird candidate crosswalk; descriptive only.
Name + region + area is a candidate, not an independently certified polygon join.
"""
import argparse, csv, io, json, math, pathlib, hashlib, re, tempfile, unicodedata, urllib.parse, zipfile
from collections import defaultdict
import numpy as np
import shapefile
from pyproj import CRS, Transformer
from tropical_pollinator_source_preflight import get, BASE, blobsha
from tropical_pollinator_source_intake import MATRICES, matrix, get_diet, norm_name, write_csv

BEE_MD5="372c22e2c24e203600a29d7f874a907d"
MODEL_SHA="a7ddae4cdb41453313d20811371cf64bd38b42b788fc89ae96c284d8d667ddd7"
SPECIES_SHA="95330e8408cf3f3816f65789c8bb374633b64342657e10051ea7e053818ac8e8"
BOX={
 "Hawaii":(-162,18,-153,24), "Lesser Antilles":(-65,10,-58,20),
 "Marianas":(143,12,148,22), "Cook Islands":(-169,-24,-155,-7),
 "Society Islands":(-156,-20,-146,-14), "Wakatobi":(122,-7.5,125,-4),
 "Cape Verde":(-26.5,14,-21.5,18), "Galapagos":(-93,-2.5,-88,2.5),
 "Riau Lingga":(102,-2.5,106,2.5), "West Sumatra":(94,-5.5,102,3),
 "Ryukyu sensitivity":(121,23,133,33),
}
ALIASES={
 ("Lesser Antilles","St.Martin"):["Saint Martin","Sint Maarten","St Martin"],
 ("Lesser Antilles","St.Barthelemy"):["Saint Barthelemy","St Barthelemy"],
 ("Lesser Antilles","St.Eustacius"):["Sint Eustatius","Saint Eustatius","St Eustatius"],
 ("Lesser Antilles","Monserrat"):["Montserrat"],
 ("Lesser Antilles","St.Lucia"):["Saint Lucia","St Lucia"],
 ("Cape Verde","Santiago2"):["Santiago"],
 ("Cape Verde","San.Antao"):["Santo Antao","Sao Antao"],
 ("Cape Verde","S.Vicente"):["Sao Vicente"],
 ("Cape Verde","S.Luzia"):["Santa Luzia"],
 ("Cape Verde","S.Nicolau"):["Sao Nicolau"],
 ("Cape Verde","Boavista"):["Boa Vista"],
 ("Galapagos","Santiago1"):["Santiago"],
 ("Galapagos","Espanola"):["Espanola"],
 ("Marianas","Aguiguan"):["Aguijan"],
 ("Marianas","Agrihan"):["Agrihan","Agrigan"],
 ("Marianas","Medinilla"):["Farallon de Medinilla"],
 ("Wakatobi","Wangi-wangi"):["Wangi Wangi","Wangi-wangi"],
 ("West Sumatra","North_Pagai"):["North Pagai","Pagai Utara"],
 ("West Sumatra","South_Pagai"):["South Pagai","Pagai Selatan"],
 ("West Sumatra","Tana_Bala"):["Tanah Bala"],
 ("West Sumatra","Tana_Massa"):["Tanah Masa"],
}
def norm(s):
    s="".join(c for c in unicodedata.normalize("NFKD",s.casefold()) if not unicodedata.combining(c))
    s=re.sub(r"\b(?:islands?|islas?|isles?|pulau)\b"," ",s)
    return "".join(c for c in s if c.isalnum())
def in_box(bounds,box):
    return all(math.isfinite(x) for x in bounds) and bounds[0]>=box[0] and bounds[1]>=box[1] and bounds[2]<=box[2] and bounds[3]<=box[3]
def load_bees():
    meta=json.loads(get("https://zenodo.org/api/records/23222636"))
    f=next(f for f in meta["files"] if f["key"]=="Data.zip")
    assert f["size"]==262554015 and f["checksum"]=="md5:"+BEE_MD5
    raw=get(f["links"].get("content") or f["links"]["self"],maxbytes=300_000_000)
    assert len(raw)==f["size"] and hashlib.md5(raw).hexdigest()==BEE_MD5
    with zipfile.ZipFile(io.BytesIO(raw)) as z, tempfile.TemporaryDirectory() as td:
        model=z.read("Data/model_table.csv")
        occurrences=z.read("Data/islandspecieslist.csv")
        assert hashlib.sha256(model).hexdigest()==MODEL_SHA
        assert hashlib.sha256(occurrences).hexdigest()==SPECIES_SHA
        rows=list(csv.DictReader(io.StringIO(model.decode("utf-8-sig"))))
        occ=list(csv.DictReader(io.StringIO(occurrences.decode("utf-8-sig"))))
        ids={r["OBJECTID_1"] for r in rows}
        assert len(rows)==len(ids)==305
        native=defaultdict(set); native_certain=defaultdict(set)
        for r in occ:
            if r["NonNative"].strip()=="0":
                native[r["OBJECTID_1"]].add(r["Corrected"])
                if r["Uncertain"].strip()=="0":
                    native_certain[r["OBJECTID_1"]].add(r["Corrected"])
        for ext in ("dbf","prj","shp","shx"):
            z.extract("Data/islands."+ext,td)
        base=pathlib.Path(td)/"Data/islands"
        crs=CRS.from_wkt(base.with_suffix(".prj").read_text())
        transform=Transformer.from_crs(crs,"EPSG:4326",always_xy=True)
        shapes=shapefile.Reader(str(base),encodingErrors="replace")
        fields=[f[0] for f in shapes.fields[1:]]
        key=next((f for f in fields if f.casefold()=="objectid_1"),None)
        if not key:
            raise RuntimeError("No exact OBJECTID_1 in source geometry: "+str(fields))
        geom={}
        for sr in shapes.iterShapeRecords():
            keyval=str(sr.record.as_dict()[key])
            if keyval.endswith(".0"): keyval=keyval[:-2]
            b=transform.transform_bounds(*sr.shape.bbox,densify_pts=21)
            if keyval in geom: raise ValueError("duplicate bee geometry ID")
            geom[keyval]={"bbox":b,"parts":len(sr.shape.parts)}
        assert ids<=geom.keys(), "missing model island geometry"
        for r in rows:
            k=r["OBJECTID_1"]
            r.update(geom[k])
        summary={"model_islands":len(rows),"geometry_islands":len(geom),"record_rows":len(occ),
                 "native_species_union":len(set().union(*native.values())),
                 "native_certain_species_union":len(set().union(*native_certain.values())),
                 "positive_model_bees":sum(float(r["FinalRichNat"])>0 for r in rows),
                 "model_native_pairs":sum(float(r["FinalRichNat"]) for r in rows),
                 "source_native_pair_count":sum(len(v) for v in native.values()),
                 "source_native_certain_pairs":sum(len(v) for v in native_certain.values()),
                 "model_richness_matches_native":sum(float(r["FinalRichNat"])==len(native[r["OBJECTID_1"]]) for r in rows),
                 "model_richness_matches_native_certain":sum(float(r["FinalRichNat"])==len(native_certain[r["OBJECTID_1"]]) for r in rows),
                 "geometry_fields":fields,"crs":crs.to_string()}
    return rows,summary

def bird_rows():
    diet,meta=get_diet()
    rows=[]
    for group,(path,expected) in MATRICES.items():
        raw=get(BASE+urllib.parse.quote(path))
        assert blobsha(raw)==expected
        names,species,area,footer=matrix(raw)
        for j,name in enumerate(names):
            ds=[diet.get(norm_name(s)) for s,v in species.items() if v[j]]
            assert all(x is not None for x in ds)
            rows.append({"group":group,"bird_name":name,"bird_area_km2":area[j]/100,
                         "all_birds":len(ds),"nectar_gt0":sum(x>0 for x in ds),
                         "nectar_ge20":sum(x>=20 for x in ds),"nectar_ge50":sum(x>=50 for x in ds)})
    return rows

def crosswalk(birds,bees):
    audit=[]
    matches=[]
    for b in birds:
        group=b["group"]; name=b["bird_name"]
        exact={norm(name)}
        alternate={norm(x) for x in ALIASES.get((group,name),[])}
        cand=[]
        if "/" not in name:
            for r in bees:
                if norm(r["IslandName"]) not in exact|alternate: continue
                if not in_box(r["bbox"],BOX[group]): continue
                ar=float(r["AreaModel"])/b["bird_area_km2"]
                cand.append((r,ar,"normalized_name" if norm(r["IslandName"]) in exact else "manual_alias"))
        eligible=[(r,ar,how) for r,ar,how in cand if .8<=ar<=1.25]
        audit.append(dict(b,status="candidate" if len(eligible)==1 else "ambiguous" if len(eligible)>1 else "unresolved",
                          name_region_candidates=len(cand),eligible_candidates=len(eligible),
                          source_candidates=";".join(r["IslandName"]+":"+r["OBJECTID_1"]+":"+str(round(ar,3)) for r,ar,_ in cand)))
        if len(eligible)!=1: continue
        r,ar,how=eligible[0]
        matches.append(dict(b,bee_id=r["OBJECTID_1"],bee_name=r["IslandName"],bee_archipelago=r["Archipelago"],
                            area_ratio=ar,match_type=how,bee_area_km2=float(r["AreaModel"]),native_bees=float(r["FinalRichNat"]),
                            dist_km=float(r["DistModel"]),dist_gt1=float(r["Dist_to_gt1"]),dist_gt1k=float(r["Dist_to_gt1k"]),
                            cov=r["Cov"],geo_origin=r["TypeModel"],biome=r["Biome"],bee_geom_parts=r["parts"],
                            bbox=";".join(str(round(x,4)) for x in r["bbox"])))
    # An island may never be counted twice, even if two independent source columns match.
    counts=defaultdict(int)
    for r in matches: counts[r["bee_id"]]+=1
    matches=[r for r in matches if counts[r["bee_id"]]==1]
    return matches,audit

def fit(rows,outcome):
    groups=sorted({r["group"] for r in rows})
    # Exploratory within-archipelago descriptive slope; no CI, no significance claims.
    if len(groups)<2 or len(rows)<len(groups)+5: return None
    x=np.array([[1.,math.log2(1+r["dist_km"]),math.log2(r["bee_area_km2"])]+[float(r["group"]==g) for g in groups[1:]] for r in rows])
    y=np.array([math.log1p(r[outcome])-math.log1p(r["native_bees"]) for r in rows])
    if np.linalg.matrix_rank(x)<x.shape[1]: return None
    beta=np.linalg.lstsq(x,y,rcond=None)[0]
    return float(beta[1])
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);p.add_argument("--self-test",action="store_true")
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    if a.self_test:
        assert norm("SÃO Vicente Island")==norm("Sao Vicente")
        assert in_box((1,1,2,2),(0,0,3,3))
        assert not in_box((1,1,4,2),(0,0,3,3))
        return
    bees,summary=load_bees()
    print("BEE_NATIVE_RECONCILIATION",json.dumps(summary),flush=True)
    birds=bird_rows()
    matches,audit=crosswalk(birds,bees)
    write_csv(a.out/"crosswalk_candidate_audit.csv",audit)
    write_csv(a.out/"candidate_bee_bird_pairs.csv",matches)
    tropical=[r for r in matches if r["group"]!="Ryukyu sensitivity" and r["native_bees"]>0]
    by_group={g:sum(r["group"]==g for r in tropical) for g in sorted({r["group"] for r in tropical})}
    descript={}
    for endpoint in ("nectar_ge50","nectar_ge20","nectar_gt0"):
        by_drop={g:fit([r for r in tropical if r["group"]!=g],endpoint) for g in by_group}
        descript[endpoint]={"within_archipelago_log1p_count_ratio_slope":fit(tropical,endpoint),"leave_one_group":by_drop}
    result={"source":summary,"candidate_pairs_all":len(matches),"candidate_pairs_tropical":len(tropical),"groups":by_group,
            "candidate_detail":[{k:r[k] for k in ("group","bird_name","bee_name","bee_id","native_bees","nectar_ge50","nectar_ge20","nectar_gt0","dist_km","area_ratio","geo_origin")} for r in matches],
            "descriptive_only":descript,"minimum_primary_support_met_on_candidate_counts":len(tropical)>=30 and len(by_group)>=10,
            "warning":"Provisional name/region/area matches, no independent bird geometry confirmation; bee-positive conditional frame; no ecological inference or confirmed cross-guild dispersal effect."}
    (a.out/"candidate_summary.json").write_text(json.dumps(result,indent=2))
    print("CANDIDATE_CROSSWALK",json.dumps(result),flush=True)
if __name__=="__main__":
    main()
