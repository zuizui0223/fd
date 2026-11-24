#!/usr/bin/env python3
"""Positive-only native-bee polygon / tropical Asian butterfly point overlap.
No inferred absence, completeness, or dispersal effect from point records.
"""
import argparse, collections, hashlib, io, json, math, pathlib, tempfile, urllib.request, zipfile
import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import shapely
from pyproj import CRS, Transformer
from tropical_pollinator_source_preflight import get
from tropical_asian_butterfly_intake import FILE_ID, MD5, BYTES
from tropical_bee_bird_crosswalk import decode_source
BEE_MD5="372c22e2c24e203600a29d7f874a907d"
ASIA=(60.,-25.,161.,30.)

def fetch_butterflies(path):
    md5=hashlib.md5();sha=hashlib.sha256();n=0
    req=urllib.request.Request("https://ndownloader.figshare.com/files/"+str(FILE_ID),headers={"User-Agent":"fd-research/1.0"})
    with urllib.request.urlopen(req,timeout=120) as r,path.open("wb") as w:
        for chunk in iter(lambda:r.read(1024*1024),b""):
            n+=len(chunk)
            if n>BYTES: raise ValueError("file size")
            md5.update(chunk);sha.update(chunk);w.write(chunk)
    assert n==BYTES and md5.hexdigest()==MD5
    return sha.hexdigest()
def fetch_islands(tmp):
    f=next(f for f in json.loads(get("https://zenodo.org/api/records/23222636"))["files"] if f["key"]=="Data.zip")
    assert f["checksum"]=="md5:"+BEE_MD5 and f["size"]==262554015
    raw=get(f["links"].get("content") or f["links"]["self"],maxbytes=300_000_000)
    assert hashlib.md5(raw).hexdigest()==BEE_MD5
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        model=pd.read_csv(io.StringIO(decode_source(z.read("Data/model_table.csv"))),dtype={"OBJECTID_1":str})
        for ext in ("dbf","prj","shp","shx"): z.extract("Data/islands."+ext,tmp)
    assert len(model)==model["OBJECTID_1"].nunique()==305
    shp=str(tmp/"Data/islands.shp")
    info=pyogrio.read_info(shp)
    crs=CRS.from_user_input(info["crs"])
    bbox=Transformer.from_crs("EPSG:4326",crs,always_xy=True).transform_bounds(*ASIA,densify_pts=21)
    data=gpd.read_file(shp,bbox=bbox,engine="pyogrio",encoding="CP1252").to_crs("EPSG:4326")
    key=next(k for k in data.columns if k.lower()=="objectid_1")
    data["bee_id"]=data[key].astype(str).str.replace(r"\.0$","",regex=True)
    assert data["bee_id"].is_unique and set(data["bee_id"])<=set(model["OBJECTID_1"])
    # Quarantine geometry crossing the broad target box, and invalid geometry.
    bounds=data.geometry.bounds.to_numpy()
    valid=data.geometry.is_valid.to_numpy() & ~data.geometry.is_empty.to_numpy()
    inside=(bounds[:,0]>=ASIA[0])&(bounds[:,1]>=ASIA[1])&(bounds[:,2]<=ASIA[2])&(bounds[:,3]<=ASIA[3])
    report={"source_crs":str(crs),"bbox_selected":len(data),"invalid_geometry":int((~valid).sum()),"not_inside_box":int((~inside).sum())}
    data=data[valid&inside].reset_index(drop=True)
    return data,model.set_index("OBJECTID_1"),report
def known(x):
    return x.notna() & ~x.str.strip().str.lower().isin(["","na","nan","none","null","unknown"])
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    with tempfile.TemporaryDirectory() as td:
        tmp=pathlib.Path(td);path=tmp/"butterfly.csv"
        sha=fetch_butterflies(path)
        geo,model,report=fetch_islands(tmp)
        print("ASIAN_BEE_GEOMETRY",json.dumps(report|{"eligible":len(geo)}),flush=True)
        if not len(geo): raise RuntimeError("no eligible Asian islands")
        tree=shapely.STRtree(geo.geometry.to_numpy())
        n=len(geo);counts=np.zeros((n,4),int);species=[[set() for _ in range(4)] for i in range(n)]
        cats={k:collections.Counter() for k in ("locality_type","how_coordinates_were_chosen","basisOfRecord")}
        total=valid_n=quality_n=ambig_n=0;all_names=set();source_valid_names=set()
        cols=["recordIndex","lohman_final_genus_species","COL_genus_species","decimalLat","decimalLon","coordinateUncertaintyInMeters","gbifID","locality_type","how_coordinates_were_chosen","basisOfRecord","reference"]
        for chunk in pd.read_csv(path,dtype=str,usecols=cols,chunksize=50000,encoding="utf-8-sig",encoding_errors="replace"):
            total+=len(chunk)
            for k in cats: cats[k].update(chunk[k].fillna("<missing>").value_counts().to_dict())
            lat=pd.to_numeric(chunk["decimalLat"],errors="coerce")
            lon=pd.to_numeric(chunk["decimalLon"],errors="coerce")
            u=pd.to_numeric(chunk["coordinateUncertaintyInMeters"],errors="coerce")
            name=chunk["lohman_final_genus_species"].str.strip()
            all_names.update(name[known(name)].unique())
            valid=lat.between(ASIA[1],ASIA[3])&lon.between(ASIA[0],ASIA[2])&known(name)
            valid_n+=int(valid.sum())
            source_valid_names.update(name[valid].unique())
            quality=valid&u.between(0,1000);quality_n+=int(quality.sum())
            frame=chunk.loc[valid].reset_index(drop=True)
            precision=quality[valid].to_numpy()
            no_gbif=(~known(frame["gbifID"])).to_numpy()
            has_reference=known(frame["reference"]).to_numpy()
            points=shapely.points(lon[valid].to_numpy(),lat[valid].to_numpy())
            q=tree.query(points,predicate="within")
            if not q.size: continue
            source,idx=q
            multiplicity=np.bincount(source,minlength=len(points))
            ambig_n+=int((multiplicity>1).sum())
            single=multiplicity[source]==1
            source,idx=source[single],idx[single]
            names=frame["lohman_final_genus_species"].str.strip().to_numpy()
            for mask_index,mask in enumerate([np.ones(len(points),bool),precision,no_gbif&has_reference,precision&no_gbif&has_reference]):
                for j in np.unique(idx[mask[source]]):
                    use=source[(idx==j)&mask[source]]
                    counts[j,mask_index]+=len(use)
                    species[j][mask_index].update(names[use])
        assert total==747825
        result=[]
        for i,row in geo.iterrows():
            r=model.loc[row["bee_id"]]
            rec={"bee_id":row["bee_id"],"island":r["IslandName"],"archipelago":r["Archipelago"],"native_bees":int(r["FinalRichNat"]),"area_km2":float(r["AreaModel"]),"distance_km":float(r["DistModel"]),"geology":r["TypeModel"]}
            for k,j in zip(("all","precision_le1km","no_gbifid_with_reference","both"),range(4)):
                rec[k+"_records"]=int(counts[i,j]);rec[k+"_recorded_species"]=len(species[i][j]) if counts[i,j]>0 else None
            result.append(rec)
        pd.DataFrame(result).to_csv(a.out/"butterfly_bee_positive_opportunity.csv",index=False)
        summary={"sha256":sha,"rows":total,"unique_final_name_strings":len(all_names),"valid_named_asia_rows":valid_n,"valid_named_asia_species_strings":len(source_valid_names),"known_uncertainty_le1km_rows":quality_n,"ambiguous_polygon_rows":ambig_n,"island_geometry":report,"eligible_bee_islands":n,"positive_island_counts":{k:sum(x[k+"_records"]>0 for x in result) for k in ("all","precision_le1km","no_gbifid_with_reference","both")},"source_category_counts":{k:v.most_common(20) for k,v in cats.items()},"warning":"Positive records only; <=1km nominal precision does not ensure uncertainty disk is inside polygon. Non-GBIF ID flag is not proof of independent source. No zeros or differential-distribution estimates."}
        (a.out/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False))
        print("ASIAN_BUTTERFLY_POLYGON_SUMMARY",json.dumps(summary,ensure_ascii=False),flush=True)
        print("ASIAN_BUTTERFLY_ISLANDS",json.dumps(result,ensure_ascii=False)[:40000],flush=True)
if __name__=="__main__": main()
