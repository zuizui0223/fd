#!/usr/bin/env python3
"""External-only EuPPollNet: same-plant susceptibility to hypothetical bee loss.
No project biological records or cohort labels; source-pinned network + geography.
"""
import argparse, hashlib, json, math, pathlib, tempfile, urllib.request, zipfile
import numpy as np
import pandas as pd

BASE="https://raw.githubusercontent.com/joseblanuza/EuPPollNet/v1.3.0/Data/3_Final_data/"
SOURCES={
 "Interaction_data.csv.gz":("cde17692dc66946002e19026433d9db429357d28873a369d884941f71131db61","6e639d4af81feee76eefa6804622f8b901c39feb"),
 "Metadata.csv":(None,"f9708e2f2b142ab2437a8be9d16e33aaefdef0ed")}
GEO="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
GEO_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
BEE={"Apidae","Andrenidae","Colletidae","Halictidae","Megachilidae","Melittidae","Stenotritidae"}
ORDERS=("hymenoptera","diptera","lepidoptera","coleoptera")
KEY=["Study_id","year","Network_id","Latitude","Longitude"]
PLANT="Plant_accepted_name"; VISITOR="Pollinator_accepted_name"
REQ=KEY[:1]+KEY[2:]+["Date","Interaction",PLANT,VISITOR,"Plant_rank","Pollinator_rank","Pollinator_order","Pollinator_family"]

def download(url,path,sha,limit,blob=None):
    h=hashlib.sha256(); n=0
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"FD source verification"}),timeout=120) as r,path.open("wb") as f:
        for b in iter(lambda:r.read(1024*1024),b""):
            n+=len(b)
            if n>limit:raise ValueError("source exceeds declared size limit")
            h.update(b);f.write(b)
    if not sha and not blob:raise ValueError("source identity pin missing")
    if sha and h.hexdigest()!=sha:raise ValueError("source SHA256 mismatch: "+path.name)
    if blob:
        raw=path.read_bytes()
        if hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()!=blob:
            raise ValueError("publisher git blob mismatch")
    return {"sha256":h.hexdigest(),"expected_sha256":sha,"bytes":n,"git_blob":blob}

def probability(n,m,d):
    if not 0<=m<=n or not 1<=d<=n:raise ValueError("invalid degree")
    return math.comb(m,d)/math.comb(n,d) if d<=m else 0.

def network_rows(g):
    pairs=g[[PLANT,VISITOR,"is_bee"]].drop_duplicates()
    if pairs.groupby(VISITOR).is_bee.nunique().gt(1).any():return None
    bees=pairs.drop_duplicates(VISITOR).set_index(VISITOR).is_bee
    n=len(bees);m=int(bees.sum())
    if n<5 or pairs[PLANT].nunique()<5:return None
    p=pairs.groupby(PLANT).is_bee.agg(["size","sum"])
    p["loss"]=(p["size"]==p["sum"]).astype(float)
    p["null"]=[probability(n,m,int(d)) for d in p["size"]]
    p["residual"]=p.loss-p["null"]
    return p.reset_index()[[PLANT,"loss","null","residual"]],n,m

def classify_areas(land):
    projected=land.to_crs(6933)
    assert projected.crs.is_projected and projected.crs.to_epsg()==6933
    assert all(a.unit_name.lower() in {"metre","meter"} for a in projected.crs.axis_info)
    area=projected.geometry.area.to_numpy()/1e6
    assert np.isfinite(area).all() and (area>0).all()
    # Never use GeoDataFrame.area attribute: it returns geometry area in
    # the current CRS even when a column named "area" has been assigned.
    return area,np.select([area>7_000_000,area>=20],["continent","island"],default="small_land")

def geography(sites,work):
    import geopandas as gpd
    z=work/"gshhg.zip";receipt=download(GEO,z,GEO_SHA,500_000_000)
    with zipfile.ZipFile(z) as a:
        members=[v for v in a.namelist() if pathlib.PurePosixPath(v).stem=="GSHHS_h_L1" and pathlib.PurePosixPath(v).suffix in {".shp",".shx",".dbf",".prj",".cpg"}]
        if len([v for v in members if v.endswith(".shp")])!=1:raise ValueError("GSHHG L1 schema drift")
        for v in members:(work/pathlib.PurePosixPath(v).name).write_bytes(a.read(v))
    # Outcome-blind known geographic controls; never enter any network.
    controls=pd.DataFrame([
      {"Study_id":"__geo_control__","year":0,"Network_id":"central_france","Latitude":46.5,"Longitude":2.5,"expected":"continent"},
      {"Study_id":"__geo_control__","year":0,"Network_id":"central_britain","Latitude":52.5,"Longitude":-1.5,"expected":"island"},
      {"Study_id":"__geo_control__","year":0,"Network_id":"central_iceland","Latitude":64.8,"Longitude":-18.5,"expected":"island"},
      {"Study_id":"__geo_control__","year":0,"Network_id":"open_atlantic","Latitude":45.0,"Longitude":-35.0,"expected":"unresolved"}])
    assert not sites.Study_id.eq("__geo_control__").any()
    all_sites=pd.concat([sites,controls[KEY]],ignore_index=True)
    box=(float(all_sites.Longitude.min())-0.1,float(all_sites.Latitude.min())-0.1,float(all_sites.Longitude.max())+0.1,float(all_sites.Latitude.max())+0.1)
    land=gpd.read_file(work/"GSHHS_h_L1.shp",bbox=box)
    if land.crs is None:land=land.set_crs(4326)
    land=land.to_crs(4326)
    idcol=next((c for c in land if c.lower()=="id"),None)
    if idcol is None:raise ValueError("GSHHG polygon IDs absent")
    land=land[[idcol,"geometry"]].dissolve(by=idcol,as_index=False)
    land["area_km2"],land["class"]=classify_areas(land)
    land["landmass"]=land[idcol].astype(str)
    points=gpd.GeoDataFrame(all_sites,geometry=gpd.points_from_xy(all_sites.Longitude,all_sites.Latitude),crs=4326)
    join=gpd.sjoin(points,land[["landmass","class","area_km2","geometry"]],how="left",predicate="within")
    bad=set(join.index[join.index.duplicated(keep=False)])
    join=join.loc[~join.index.isin(bad),KEY+["landmass","class","area_km2"]].copy()
    join["class"]=join["class"].fillna("unresolved")
    out=all_sites.merge(join,on=KEY,how="left",validate="one_to_one")
    out["class"]=out["class"].fillna("unresolved")
    check=controls.merge(out[["Network_id","class","area_km2"]].loc[out.Study_id.eq("__geo_control__")],on="Network_id",validate="one_to_one")
    receipt["fixed_geography_controls"]=[{"control":v.Network_id,"expected":v.expected,"observed":v["class"],"area_km2":float(v.area_km2) if pd.notna(v.area_km2) else None} for _,v in check.iterrows()]
    print("FD_EXTERNAL_GEO_CONTROLS",json.dumps(receipt["fixed_geography_controls"],allow_nan=False),flush=True)
    assert len(check)==len(controls) and check.expected.eq(check["class"]).all(),"positive geographic controls failed"
    receipt["area_crs"]="EPSG:6933; projected square metres divided by 1e6"
    receipt["landmass_classes"]={str(k):int(v) for k,v in land["class"].value_counts().items()}
    out=out.loc[~out.Study_id.eq("__geo_control__")].drop(columns="area_km2").copy()
    out["block"]=[f"{math.floor(y/5)}:{math.floor(x/5)}" for y,x in zip(out.Latitude,out.Longitude)]
    receipt["ambiguous_sites_excluded"]=len(bad)
    return out,receipt

def analyze(d,meta,work):
    assert set(REQ).issubset(d),sorted(set(REQ)-set(d))
    assert {"Study_id","Year","Taxa_recorded"}.issubset(meta)
    census={"source_rows":len(d),"source_studies":int(d.Study_id.nunique()),
       "plant_ranks":{str(k):int(v) for k,v in d.Plant_rank.fillna("[missing]").value_counts().items()},
       "visitor_ranks":{str(k):int(v) for k,v in d.Pollinator_rank.fillna("[missing]").value_counts().items()}}
    bad=d[REQ].fillna("").apply(lambda x:x.str.contains("\ufffd",regex=False)).any(axis=1)
    if bad.mean()>.01:raise ValueError("invalid source encoding exceeds 1%")
    d=d.loc[~bad,REQ].copy()
    d["year"]=pd.to_numeric(d.Date.str.extract(r"^(\d{4})-",expand=False),errors="coerce")
    meta=meta.copy();meta["year"]=pd.to_numeric(meta.Year,errors="coerce")
    meta["broad"]=meta.Taxa_recorded.fillna("").str.lower().map(lambda v:all(k in v for k in ORDERS))
    broad=meta.groupby(["Study_id","year"]).broad.agg(["all","any"]).reset_index()
    broad=broad.loc[broad["all"]&broad["any"],["Study_id","year"]]
    d=d.merge(broad,on=["Study_id","year"],how="inner",validate="many_to_one")
    census["four_order_documented_rows"]=len(d)
    for c in ("Latitude","Longitude"):d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.loc[d.Latitude.between(-90,90)&d.Longitude.between(-180,180)&
      d.Plant_rank.str.lower().eq("species")&d.Pollinator_rank.str.lower().eq("species")&
      d[PLANT].notna()&d[VISITOR].notna()&d.Pollinator_order.notna()&d.Pollinator_family.notna()&
      pd.to_numeric(d.Interaction,errors="coerce").gt(0)&d.Network_id.notna()].copy()
    d["is_bee"]=d.Pollinator_order.eq("Hymenoptera")&d.Pollinator_family.isin(BEE)
    census["qualified_visit_rows"]=len(d);census["qualified_studies"]=int(d.Study_id.nunique())
    if d.empty:return {"census":census,"status":"no_qualified_records"}
    sites=d[KEY].drop_duplicates().reset_index(drop=True)
    loc,receipt=geography(sites,work)
    census["site_classes"]={str(k):int(v) for k,v in loc["class"].value_counts().items()}
    # Same year and exact coordinate, not an invented temporal replicate.
    x=d.merge(loc,on=KEY,how="inner",validate="many_to_one")
    rec=[];net=[]
    for key,g in x.groupby(KEY,sort=True):
        if g["class"].iloc[0] not in ("island","continent"):continue
        q=network_rows(g)
        if q is None:continue
        p,n,m=q
        for k,v in zip(KEY,key):p[k]=v
        for k in ("class","landmass","block"):p[k]=g[k].iloc[0]
        rec.append(p)
        net.append({"class":g["class"].iloc[0],"study":key[0],"n":n,"m":m,"res":float(p.residual.mean())})
    if not rec:return {"census":census,"geography":receipt,"status":"no_eligible_networks"}
    p=pd.concat(rec,ignore_index=True);net=pd.DataFrame(net)
    # Equal network-site weights within a plant x study x year x land class.
    a=p.groupby(["Study_id","year",PLANT,"class"])[["loss","null","residual"]].mean().unstack("class")
    needed=[(k,c) for k in ("loss","null","residual") for c in ("island","continent")]
    if not all(k in a for k in needed):a=pd.DataFrame(columns=pd.MultiIndex.from_tuples(needed))
    else:a=a.dropna(subset=needed)
    matched=p.merge(a.index.to_frame(index=False),on=["Study_id","year",PLANT],how="inner",validate="many_to_one") if not a.empty else p.iloc[:0]
    support={"matched_plant_study_year":len(a),"matched_species":int(matched[PLANT].nunique()),"matched_studies":int(matched.Study_id.nunique()),
      "matched_island_landmasses":int(matched.loc[matched["class"].eq("island"),"landmass"].nunique()),
      "matched_blocks":{c:int(matched.loc[matched["class"].eq(c),"block"].nunique()) for c in ("island","continent")}}
    gate=support["matched_studies"]>=8 and support["matched_species"]>=20 and support["matched_island_landmasses"]>=5 and min(support["matched_blocks"].values())>=5
    result={"census":census,"geography":receipt,"network_support":{c:{"networks":int(len(g)),"studies":int(g.study.nunique()),"zero_bee":int(g.m.eq(0).sum()),"all_bee":int(g.m.eq(g.n).sum()),"descriptive_network_mean_residual":float(g.res.mean())} for c,g in net.groupby("class")},
      "same_plant_support":support,"primary_support_gate":bool(gate),"primary_supported":False}
    if not a.empty:
        delta=pd.DataFrame({k:a[(k,"island")]-a[(k,"continent")] for k in ("loss","null","residual")})
        by=delta.groupby(level="Study_id").mean()
        result["descriptive_matched_island_minus_continent"]={k:float(v) for k,v in by.mean().items()}
        if gate:
            v=by.residual.to_numpy();rng=np.random.default_rng(20261128)
            boot=v[rng.integers(len(v),size=(9999,len(v)))].mean(1)
            loo=(v.sum()-v)/(len(v)-1);ci=np.quantile(boot,[.025,.975])
            result["study_bootstrap_residual_interval"]=ci.tolist()
            result["leave_study_residual_range"]=[float(loo.min()),float(loo.max())]
            result["primary_supported"]=bool((ci[0]>0 and loo.min()>0) or (ci[1]<0 and loo.max()<0))
    return result

def self_test():
    assert probability(5,2,2)==.1 and probability(5,2,3)==0
    # Exhaustive equal-size removal confirms the exact hypergeometric probability.
    import itertools
    for n in range(2,8):
        for m in range(n+1):
            for d in range(1,n+1):
                s=list(itertools.combinations(range(n),m))
                assert abs(probability(n,m,d)-sum(set(range(d)).issubset(x) for x in s)/len(s))<1e-12
    g=pd.DataFrame({PLANT:["p"+str(i) for i in range(5) for _ in range(2)],VISITOR:["v"+str(j) for i in range(5) for j in (i,(i+1)%5)]})
    g["is_bee"]=g[VISITOR].isin(["v0","v1"])
    p,n,m=network_rows(g);assert n==5 and m==2 and abs(p["null"].mean()-.1)<1e-12
    # Geometry regression: a large and a small projected polygon with an
    # intentionally conflicting "area" attribute in the table.
    import geopandas as gpd
    from shapely.geometry import box
    fake=gpd.GeoDataFrame({"area":[0.,99999999.]},geometry=[box(0,0,3000000,3000000),box(0,0,10000,10000)],crs=6933).to_crs(4326)
    ar,classes=classify_areas(fake)
    assert np.allclose(ar,[9000000,100],rtol=1e-5) and list(classes)==["continent","island"]
    # Integration fixture: same plants, two classes, only one study;
    # support must fail, but the conditional estimate must be computable.
    from unittest.mock import patch
    visits=[]
    for role,lat in (("island","31"),("continent","41")):
        for i in range(5):
            for j in (i,(i+1)%5):
                visits.append({"Study_id":"s","Network_id":role,"Latitude":lat,"Longitude":"0",
                   "Date":"2020-06-01","Interaction":"1",PLANT:"Plant "+str(i),VISITOR:"Visitor "+str(j),
                   "Plant_rank":"species","Pollinator_rank":"species",
                   "Pollinator_order":"Hymenoptera" if j<2 else "Diptera",
                   "Pollinator_family":"Apidae" if j<2 else "Syrphidae"})
    metadata=pd.DataFrame({"Study_id":["s"],"Year":["2020"],"Taxa_recorded":["Hymenoptera, Diptera, Lepidoptera and Coleoptera"]})
    def fake_geo(sites,work):
        sites=sites.copy()
        sites["class"]=sites.Network_id
        sites["landmass"]=sites.Network_id
        sites["block"]=sites.Network_id
        return sites,{"synthetic":True}
    with patch(__name__+".geography",side_effect=fake_geo):
        q=analyze(pd.DataFrame(visits),metadata,None)
    assert q["same_plant_support"]["matched_species"]==5
    assert not q["primary_support_gate"] and not q["primary_supported"]
    assert abs(q["descriptive_matched_island_minus_continent"]["residual"])<1e-12
    json.dumps(q,allow_nan=False)
    print("FD_EXTERNAL_BEE_LOSS_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true");p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    self_test()
    if a.self_test:return
    with tempfile.TemporaryDirectory() as td:
        w=pathlib.Path(td);pins={}
        for f,(sha,blob) in SOURCES.items():pins[f]=download(BASE+f,w/f,sha,5_000_000,blob)
        d=pd.read_csv(w/"Interaction_data.csv.gz",dtype=str,encoding="utf-8",encoding_errors="replace",low_memory=False)
        if len(d)!=623476:raise ValueError("pinned source row count mismatch")
        m=pd.read_csv(w/"Metadata.csv",dtype=str)
        result=analyze(d,m,w)
        result["source_pins"]=pins
        result["scope"]="External source only, observed species-level visits, documented four-order coverage not guaranteed complete survey. Hypothetical deletion and conditional random-species expectation, not measured pollination or historical loss."
        result["primary"]="Study-equal same-plant within-study/year island-minus-continent excess disconnection; support gate fixed before outcomes. Geography class is positive polygon membership, no inference from absence."
        a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
        print("FD_EXTERNAL_BEE_LOSS_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
