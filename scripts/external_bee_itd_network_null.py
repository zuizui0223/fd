#!/usr/bin/env python3
"""External-only EuPPollNet x EuPollTrait bee morphology interaction test.

No project botanical data. Whole local observed bee subwebs, fixed binary
degree sequences, equal-plant breadth. ITD is measured body size, not flight
ability, proboscis length, or pollination effectiveness.
"""
import argparse
import collections
import hashlib
import io
import json
import pathlib
import re
import subprocess
import tempfile
import zipfile

import numpy as np
import pandas as pd
import openpyxl
from numba import njit

from tropical_pollinator_source_preflight import get

EU_REF="54bcb8e23bdc221d930a6d212778103bbbff5287"
EU_BASE="https://raw.githubusercontent.com/joseblanuza/EuPPollNet/"+EU_REF+"/Data/3_Final_data/"
EU_FILES={"Interaction_data.csv.gz":"6e639d4af81feee76eefa6804622f8b901c39feb","Metadata.csv":"f9708e2f2b142ab2437a8be9d16e33aaefdef0ed"}
TRAIT_URL="https://zenodo.org/records/18032357/files/pollinator_traits_bees_hoverflies_repository.zip?download=1"
TRAIT_MD5="7a50cc03e2c5a287e663dc5d029f8026"
GEO_URL="https://github.com/GenericMappingTools/gshhg-gmt/releases/download/2.3.7/gshhg-shp-2.3.7.zip"
GEO_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
SEED=20261128
DRAWS=399

def blob(b):
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()

def digest(b):
    return hashlib.sha256(b).hexdigest()

def clean(x):
    if x is None:
        return ""
    return " ".join(str(x).replace("_"," ").strip().split())

def name(x):
    x=clean(x)
    if re.fullmatch(r"[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+",x):
        return x
    return None

def num(x):
    try:
        y=float(x)
        return y if np.isfinite(y) and y>0 else np.nan
    except (TypeError,ValueError):
        return np.nan

def texttable(raw,compression=None):
    try:
        return pd.read_csv(io.BytesIO(raw),compression=compression,dtype=str,keep_default_na=False,encoding="utf-8-sig")
    except UnicodeDecodeError:
        return pd.read_csv(io.BytesIO(raw),compression=compression,dtype=str,keep_default_na=False,encoding="cp1252")

def traits():
    raw=get(TRAIT_URL,maxbytes=10_000_000)
    assert hashlib.md5(raw).hexdigest()==TRAIT_MD5,"EuPollTrait publisher checksum mismatch"
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        fn=[x for x in z.namelist() if x.endswith("/output/traits_raw_bees_hoverflies.xlsx") or x=="output/traits_raw_bees_hoverflies.xlsx"]
        if len(fn)!=1:
            raise ValueError("expected one published raw workbook: "+str(fn))
        md=[x for x in z.namelist() if x.endswith("/output/metadata_traits_raw.xlsx") or x=="output/metadata_traits_raw.xlsx"]
        assert len(md)==1,"publisher metadata workbook absent"
        wb=openpyxl.load_workbook(io.BytesIO(z.read(fn[0])),read_only=True,data_only=True)
        mb=openpyxl.load_workbook(io.BytesIO(z.read(md[0])),read_only=True,data_only=True)
        ws=wb.active
        rows=list(ws.values)
        meta=[list(map(clean,row)) for sheet in mb.worksheets for row in sheet.values]
        wb.close();mb.close()
    assert len(rows)>=3000 and len(rows[0])>=100,"unexpected workbook dimensions"
    # Publisher data descriptor defines first five columns as order, family,
    # genus, epithet, combined name; first three rows are hierarchical headers.
    headers=[list(map(clean,row)) for row in rows[:3]]
    parents=[[],[]]
    for level in range(2):
        last=""
        for s in headers[level]:
            if s:
                last=s
            parents[level].append(last)
    candidates=[]
    for j in range(5,len(headers[2])):
        label=" | ".join([parents[0][j],parents[1][j],headers[2][j]])
        # Use a specific female/worker value, not male/queen, N, range or
        # derived tongue length. Never auto-select among ambiguous columns.
        lo=label.casefold()
        leaf=headers[2][j].casefold()
        if (re.search(r"\bitd\b|inter[ -]?tegular",lo)
            and not re.search(r"tongue|hair",lo)
            and re.search(r"\bfemale\b|\bworker\b|♀",lo)
            and not re.search(r"\bmale\b|\bqueen\b|♂",leaf)
            and not re.search(r"\b(n|count|number|sample|sd|se|min|max)\b|specimen|deviation",leaf)):
            candidates.append((j,label))
    relevant_meta=[r for r in meta if any(re.search(r"\bitd\b|inter[ -]?tegular",s,re.I) for s in r)]
    info={"archive_sha256":digest(raw),"publisher_md5":TRAIT_MD5,"dimensions":[len(rows),len(rows[0])],"headers":headers,"itd_metadata":relevant_meta,"female_candidates":candidates}
    print("FD_EUPOLLTRAIT_SCHEMA",json.dumps(info,ensure_ascii=True),flush=True)
    if len(candidates)!=1:
        return {},dict(info,status="ambiguous_female_itd_no_fit")
    j,label=candidates[0]
    # Require source metadata to corroborate the millimetre unit.
    if not any(re.search(r"\bmm\b|millimet",s,re.I) for r in relevant_meta for s in r):
        return {},dict(info,status="itd_unit_not_confirmed_no_fit")
    by=collections.defaultdict(set)
    for r in rows[3:]:
        if len(r)<=j or clean(r[0]).casefold()!="hymenoptera":
            continue
        n=name(r[4])
        # Require redundant genus and epithet to match exactly.
        if not n or n!=name(clean(r[2])+" "+clean(r[3])):
            continue
        v=num(r[j])
        if np.isfinite(v) and .1<=v<=15:
            by[n].add(v)
    out={k:float(np.log2(next(iter(v)))) for k,v in by.items() if len(v)==1}
    return out,dict(info,status="verified_female_itd" if out else "no_valid_values",selected_column=j,selected_label=label,valid_species=len(out),conflicting_names=sum(len(v)>1 for v in by.values()))

def source_visits():
    data={}
    pins={}
    for f,h in EU_FILES.items():
        b=get(EU_BASE+f,maxbytes=20_000_000)
        assert blob(b)==h,"EuPPollNet source blob mismatch"
        data[f]=texttable(b,"gzip" if f.endswith(".gz") else None)
        pins[f]={"blob":h,"sha256":digest(b),"rows":len(data[f])}
    d=data["Interaction_data.csv.gz"]
    assert len(d)==623476,"EuPPollNet v1.3.0 rows changed"
    req=["Study_id","Network_id","Latitude","Longitude","Plant_accepted_name","Pollinator_accepted_name","Date"]
    assert set(req)<=set(d.columns),"source schema drift"
    d=d[req].copy()
    d["plant"]=d.Plant_accepted_name.map(name)
    d["visitor"]=d.Pollinator_accepted_name.map(name)
    # Source ISO year, never guess an ambiguous date.
    d["year"]=d.Date.str.extract(r"^((?:19|20)\d{2})(?:[-/]|$)",expand=False)
    d["lat"]=pd.to_numeric(d.Latitude,errors="coerce")
    d["lon"]=pd.to_numeric(d.Longitude,errors="coerce")
    keep=d.plant.notna()&d.visitor.notna()&d.year.notna()&d.lat.between(-90,90)&d.lon.between(-180,180)&d.Study_id.ne("")&d.Network_id.ne("")
    census={"raw_rows":len(d),"qualified_name_date_coord_rows":int(keep.sum()),"bad_dates":int(d.year.isna().sum()),"source_pins":pins,"study_count":int(d.Study_id.nunique())}
    print("FD_EXTERNAL_EU_NETWORK_CENSUS",json.dumps(census),flush=True)
    d=d.loc[keep,["Study_id","Network_id","year","lat","lon","plant","visitor"]].drop_duplicates()
    return d,census

@njit
def swaps(edges,A,n,seed):
    np.random.seed(seed)
    moved=0
    for _ in range(n):
        i,j=np.random.randint(len(edges)),np.random.randint(len(edges))
        p,b=edges[i,0],edges[i,1]
        q,c=edges[j,0],edges[j,1]
        if p==q or b==c or A[p,c] or A[q,b]:
            continue
        A[p,b]=False;A[q,c]=False
        A[p,c]=True;A[q,b]=True
        edges[i,1]=c;edges[j,1]=b
        moved+=1
    return moved

@njit
def breadth(A,values):
    total=0.;n=0
    for p in range(A.shape[0]):
        ids=np.where(A[p])[0]
        k=len(ids)
        if k<2:
            continue
        s=0.
        for i in range(k):
            for j in range(i):
                s+=abs(values[ids[i]]-values[ids[j]])
        total+=s/(k*(k-1)/2)
        n+=1
    return total/n if n else np.nan

def fit_networks(d,traits):
    out=[]
    reject=collections.Counter()
    # Require one exact coordinate per original source study/site/year,
    # not multiple sites or inferred centroids.
    for k,df in d.groupby(["Study_id","Network_id","year"],sort=True):
        loc=df[["lat","lon"]].drop_duplicates()
        if len(loc)!=1:
            reject["multi_coord_site_year"]+=1
            continue
        df=df[df.visitor.isin(traits)]
        if df.empty:
            reject["no_trait_visitor"]+=1
            continue
        P=sorted(set(df.plant));B=sorted(set(df.visitor))
        if len(P)<5 or len(B)<5:
            reject["small_subweb"]+=1
            continue
        pi={x:i for i,x in enumerate(P)};bi={x:i for i,x in enumerate(B)}
        A=np.zeros((len(P),len(B)),dtype=np.bool_)
        A[[pi[x] for x in df.plant],[bi[x] for x in df.visitor]]=True
        if (A.sum(1)>=2).sum()<3 or len({traits[x] for x in B})<3:
            reject["insufficient_breadth"]+=1
            continue
        E=np.argwhere(A).astype(np.int64);n=len(E)
        vals=np.array([traits[x] for x in B],dtype=float)
        d0=A.sum(0).copy();d1=A.sum(1).copy()
        observed=breadth(A,vals)
        seed=(SEED+int(hashlib.sha256("|".join(map(str,k)).encode()).hexdigest()[:8],16))%(2**31-1)
        moved=swaps(E,A,20*n,seed)
        if moved<max(10,n):
            reject["burnin_low_switch"]+=1
            continue
        draw=[];switch=[]
        for i in range(DRAWS):
            switch.append(swaps(E,A,5*n,(seed+i+1)%(2**31-1)))
            draw.append(breadth(A,vals))
        assert np.array_equal(A.sum(0),d0) and np.array_equal(A.sum(1),d1)
        if np.mean(switch)<max(2,.1*n) or np.std(draw)<1e-12:
            reject["low_switch_or_degenerate_null"]+=1
            continue
        lat,lon=loc.iloc[0]
        out.append({"study":str(k[0]),"lat":float(lat),"lon":float(lon),"block":(int(np.floor(lat/10)),int(np.floor(lon/10))),"plants":len(P),"bees":len(B),"links":n,"observed":float(observed),"draw":np.array(draw),"mean_switches":float(np.mean(switch))})
    return out,dict(reject)

def geography(rr):
    import geopandas as gpd
    from shapely import make_valid
    from shapely.ops import unary_union
    from shapely.geometry import Point
    raw=get(GEO_URL,maxbytes=300_000_000)
    assert digest(raw)==GEO_SHA,"GSHHG pinned source checksum mismatch"
    with tempfile.TemporaryDirectory() as tmp,zipfile.ZipFile(io.BytesIO(raw)) as z:
        dest=pathlib.Path(tmp)
        found=[]
        for n in z.namelist():
            p=pathlib.PurePosixPath(n)
            if p.name.startswith("GSHHS_h_L1."):
                (dest/p.name).write_bytes(z.read(n));found.append(p.name)
        assert "GSHHS_h_L1.shp" in found
        g=gpd.read_file(dest/"GSHHS_h_L1.shp")
    if g.crs is None:
        raise ValueError("GSHHG CRS absent")
    g=g.to_crs(4326)
    key=next((x for x in ("sibling_id","SIBLING_ID","id","ID") if x in g),None)
    if key is None:
        raise ValueError("GSHHG source id absent")
    grouped=g.groupby(key,sort=True)
    pieces=[(str(i),make_valid(unary_union(list(v.geometry)))) for i,v in grouped]
    land=gpd.GeoDataFrame({"unit":[a for a,_ in pieces]},geometry=[b for _,b in pieces],crs=4326)
    land["area_km2"]=land.to_crs(6933).area/1e6
    # Positive membership only, no unmatched point treated as mainland.
    # Include islands >=20 km2; smaller/ambiguous/unmapped remain unresolved.
    land=land.loc[land.area_km2>=20].copy()
    assert 3<=int((land.area_km2>7_000_000).sum())<=9
    pts=gpd.GeoDataFrame({"ix":range(len(rr))},geometry=[Point(z["lon"],z["lat"]) for z in rr],crs=4326)
    j=gpd.sjoin(pts,land[["unit","area_km2","geometry"]],predicate="within",how="inner")
    for i,z in enumerate(rr):
        m=j.loc[j.ix.eq(i)]
        if len(m)!=1:
            z["class"]="unresolved"
            continue
        a=m.iloc[0]
        z["class"]="mainland" if a.area_km2>7_000_000 else "island"
        z["unit"]=str(a.unit)
        geom=land.loc[land.unit.eq(a.unit),"geometry"].iloc[0]
        z["coast_clear"]=bool(geom.boundary.distance(pts.geometry.iloc[i])>.005)
    return {"source_sha256":digest(raw),"area_min_km2":20,"continental_threshold_km2":7_000_000,"unmatched_or_small_is_unresolved":True,"classification_counts":dict(collections.Counter(x["class"] for x in rr))}

def grouped(rr,key):
    blocks=sorted({r[key] for r in rr})
    o=[];d=[]
    for b in blocks:
        sub=[r for r in rr if r[key]==b]
        o.append(np.mean([r["observed"] for r in sub]))
        d.append(np.mean([r["draw"] for r in sub],axis=0))
    return np.array(o),np.array(d)

def pvalue(o,draw):
    return min(1.,2*min((1+np.sum(draw<=o+1e-12))/(1+len(draw)),(1+np.sum(draw>=o-1e-12))/(1+len(draw))))

def summary(rr,key):
    if not rr:
        return {"networks":0,"groups":0,"status":"no_support"}
    o,d=grouped(rr,key)
    y=o-d.mean(1);rng=np.random.default_rng(SEED)
    lo=[np.delete(y,i).mean() for i in range(len(y))] if len(y)>1 else []
    boot=[rng.choice(y,len(y),replace=True).mean() for _ in range(999)]
    return {"networks":len(rr),"groups":len(y),"studies":len({r["study"] for r in rr}),"observed":float(o.mean()),"null_mean":float(d.mean()),"residual":float(y.mean()),"conditional_p":pvalue(o.mean(),d.mean(0)),"group_boot_interval":np.quantile(boot,[.025,.975]).tolist() if len(y)>=10 else None,"leave_group_range":[float(min(lo)),float(max(lo))] if lo else None,"support_gate":len(rr)>=20 and len(y)>=10}

def contrast(rr):
    a=[r for r in rr if r.get("class")=="island"]
    b=[r for r in rr if r.get("class")=="mainland"]
    aa=summary(a,"block");bb=summary(b,"block")
    out={"island":aa,"mainland":bb,"island_units":len({r["unit"] for r in a}),"mixed_studies":len({r["study"] for r in a}&{r["study"] for r in b})}
    gate=aa.get("support_gate",False) and bb.get("support_gate",False) and aa.get("studies",0)>=3 and bb.get("studies",0)>=3 and out["island_units"]>=10 and out["mixed_studies"]>=5
    out["support_gate"]=bool(gate)
    # Do not run regional inference below the gate. Show descriptive means.
    if a and b:
        out["descriptive_difference"]=aa["residual"]-bb["residual"]
    if gate:
        # Same geographic-block multiplier applied across both classes; do
        # not pretend overlapping geographic blocks are independent.
        keys=sorted({r["block"] for r in a+b})
        oa,da=grouped(a,"block");ob,db=grouped(b,"block")
        ka=sorted({r["block"] for r in a});kb=sorted({r["block"] for r in b})
        ya=oa-da.mean(1);yb=ob-db.mean(1)
        rng=np.random.default_rng(SEED)
        boot=[]
        for _ in range(999):
            w=dict(zip(keys,rng.exponential(size=len(keys))))
            boot.append(np.average(ya,weights=[w[k] for k in ka])-np.average(yb,weights=[w[k] for k in kb]))
        loo=[]
        for k in keys:
            ma=np.array([x!=k for x in ka]);mb=np.array([x!=k for x in kb])
            if ma.any() and mb.any():loo.append(ya[ma].mean()-yb[mb].mean())
        out.update(conditional_p=pvalue(oa.mean()-ob.mean(),da.mean(0)-db.mean(0)),joint_group_range=np.quantile(boot,[.025,.975]).tolist(),leave_group_range=[float(min(loo)),float(max(loo))])
        # Within-study geographical contrast, descriptive check of sampling
        # study confounding; no claim of historical island effect.
        dif=[]
        for s in sorted({r["study"] for r in a}&{r["study"] for r in b}):
            xa=summary([r for r in a if r["study"]==s],"block")
            xb=summary([r for r in b if r["study"]==s],"block")
            dif.append(xa["residual"]-xb["residual"])
        out["within_study_mean_difference"]=float(np.mean(dif))
    return out

def self_test():
    A=np.array([[1,1,0],[0,1,1],[1,0,1]],dtype=bool);v=np.array([0.,1.,3.])
    assert abs(breadth(A,v)-2)<1e-12
    E=np.argwhere(A);a=A.sum(0).copy();b=A.sum(1).copy()
    assert swaps(E,A,1000,SEED)>0
    assert np.array_equal(a,A.sum(0)) and np.array_equal(b,A.sum(1))
    assert name("Bombus_terrestris")=="Bombus terrestris" and name("Bombus sp.") is None
    print("FD_BEE_ITD_NULL_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=pathlib.Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args();self_test()
    if a.self_test:return
    if a.out is None:p.error("--out required")
    t,ti=traits();d,di=source_visits()
    di["matched_trait_visitor_names"]=len(set(d.visitor)&t.keys())
    di["matched_trait_rows"]=int(d.visitor.isin(t).sum())
    result={"trait_source":ti,"network_source":di,"fit_status":"not_run","geography":None}
    if t:
        rr,rej=fit_networks(d,t)
        result["fit_status"]="completed"
        result["rejections"]=rej
        result["all_block"]=summary(rr,"block")
        result["all_study"]=summary(rr,"study")
        if rr:
            try:
                result["geography"]=geography(rr)
                result["island_mainland"]=contrast(rr)
                result["coast_clear_sensitivity"]=contrast([r for r in rr if r.get("coast_clear",False)])
            except Exception as exc:
                result["geography"]={"status":"unavailable_no_regional_fit","error":str(exc)[:600]}
    # Global and geographic conditional nulls are a fixed two-endpoint family;
    # no support below gates. No result promotion if geography fails.
    p0=result.get("all_block",{}).get("conditional_p",1) if result.get("all_block",{}).get("support_gate",False) else 1.
    p1=result.get("island_mainland",{}).get("conditional_p",1) if result.get("island_mainland",{}).get("support_gate",False) else 1.
    pvals=[p0,p1];order=np.argsort(pvals);adj=[1.,1.];last=0.
    for i,j in enumerate(order):
        last=max(last,min(1.,(2-i)*pvals[j]));adj[j]=last
    result["holm2"]={"global_block":adj[0],"island_mainland":adj[1]}
    result["claim_limit"]="Exact-name recorded bee subwebs only. ITD is adult size; degree-conditioned breadth not effective pollination, nectar accessibility or historical filtering. Shared species/studies and source coverage limit conditional uncertainty; no inferred absences. Geography uses positive source polygon membership, no own biological cohort."
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False,default=lambda x:x.item() if hasattr(x,"item") else str(x)))
    print("FD_BEE_ITD_NULL_RESULT",json.dumps({k:v for k,v in result.items() if k!="trait_source"},allow_nan=False,default=str),flush=True)
if __name__=="__main__":
    main()
