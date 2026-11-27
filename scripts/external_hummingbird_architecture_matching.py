#!/usr/bin/env python3
"""Cross-network floral symmetry x tube-bill matching on actual published links.

Marjakangas et al. 2025 Dryad data reuse (79 networks from Dalsgaard 2021).
Independent NHM corolla tube measurements and Wang symmetry, no imputed lengths.
Degree-preserving randomization conditional on retained measured complete subwebs.
No inference of effective pollination or island identity from coordinates.
"""
import argparse,collections,hashlib,io,json,pathlib,re,subprocess,tempfile,urllib.parse
import numpy as np
import pandas as pd
from tropical_pollinator_source_preflight import get
from external_floral_replication_intake import name
from external_gift_floral_replication_coverage import source_tables
DOI="10.5061/dryad.h9w0vt4sq"
API="https://datadryad.org/api/v2"
SEED=20261128

def fetch_json(url):
    b=get(url,maxbytes=5_000_000)
    return json.loads(b),hashlib.sha256(b).hexdigest()

def absolute(href):
    if href.startswith("/"):url="https://datadryad.org"+href
    else:url=href
    u=urllib.parse.urlparse(url)
    assert u.scheme=="https" and u.hostname=="datadryad.org","unexpected publisher URL"
    return url

def fetch_network():
    meta,md=fetch_json(API+"/datasets/"+urllib.parse.quote("doi:"+DOI,safe=""))
    assert DOI in str(meta.get("identifier","")).lower()
    href=meta["_links"]["stash:version"]["href"]
    ver,vd=fetch_json(absolute(href))
    files,fd=fetch_json(absolute(ver["_links"]["stash:files"]["href"]))
    fs=files["_embedded"]["stash:files"]
    f=[x for x in fs if pathlib.PurePosixPath(x["path"]).name=="network_long.RData"]
    assert len(f)==1,f
    f=f[0];raw=get(absolute(f["_links"]["stash:download"]["href"]),maxbytes=5_000_000)
    assert len(raw)==int(f["size"])
    alg=str(f.get("digestType","")).lower().replace("-","");expected=str(f.get("digest","")).lower()
    assert alg in {"md5","sha256","sha1"} and hashlib.new(alg,raw).hexdigest()==expected,"publisher digest unavailable or mismatch"
    with tempfile.TemporaryDirectory() as td:
        p=pathlib.Path(td);(p/"net.RData").write_bytes(raw)
        expr='a<-commandArgs(TRUE);e<-new.env(parent=emptyenv());nm<-load(a[1],envir=e);need<-c("Plant_species","Animal_species","Network_name","Latitude","Longitude","Beak.Length_Culmen","outcome");keep<-nm[sapply(nm,function(x)is.data.frame(e[[x]])&&all(need%in%names(e[[x]])))];stopifnot(length(keep)==1);write.csv(e[[keep]],a[2],row.names=FALSE,na="")'
        subprocess.run(["Rscript","--vanilla","-e",expr,str(p/"net.RData"),str(p/"net.csv")],check=True)
        d=pd.read_csv(p/"net.csv",dtype=str,keep_default_na=False)
    info={"doi":DOI,"dataset_metadata_sha256":md,"version_metadata_sha256":vd,"file_metadata_sha256":fd,"version":ver.get("version"),"path":f["path"],"size":len(raw),"publisher_digest":expected,"digest_type":alg,"sha256":hashlib.sha256(raw).hexdigest(),"rows":len(d),"columns":list(d.columns)}
    print("FD_HUMMINGBIRD_SOURCE",json.dumps(info),flush=True)
    return d,info

def binomial(v):
    if not isinstance(v,str):return None
    return name(" ".join(v.strip().replace("_"," ").split()))

def swap(E,A,rng,attempts):
    # Select two existing links uniformly; only legal checkerboard switches.
    # Each switch preserves both degree sequences and the total edge count.
    moved=0
    for _ in range(attempts):
        a,b=rng.integers(len(E),size=2)
        i,j=E[a];u,v=E[b]
        if i==u or j==v or A[i,v] or A[u,j]:continue
        A[i,j]=A[u,v]=False;A[i,v]=A[u,j]=True
        E[a,1]=v;E[b,1]=j;moved+=1
    return moved

def metrics(A,L,radial):
    # Equal weight to plants, not links or bird abundance. Compare observed
    # mean absolute log2 bill/tube ratio within symmetry classes.
    y=(A*L).sum(1)/A.sum(1)
    return np.array([y[radial].mean(),y[~radial].mean(),y[radial].mean()-y[~radial].mean()])

def panels(d,tab,draws=399):
    out=[];reject=collections.Counter();rng=np.random.default_rng(SEED)
    for net,rows in d.groupby("Network_name",sort=True):
        if not net:reject["unnamed"]+=1;continue
        if rows.duplicated(["p","a"]).any():reject["duplicate_pair"]+=1;continue
        P=sorted(set(rows.p));B=sorted(set(rows.a))
        # Complete pair universe before conditioning on known measured traits.
        if len(rows)!=len(P)*len(B):reject["incomplete_pair_universe"]+=1;continue
        rows=rows[rows.p.isin(tab)&rows.bill.notna()].copy()
        if rows.empty:reject["no_measured_carriers"]+=1;continue
        A=rows.pivot(index="p",columns="a",values="y")
        if A.isna().any().any():reject["incomplete_measured_subweb"]+=1;continue
        A=A.loc[A.sum(axis=1)>0,A.sum(axis=0)>0]
        if A.shape[0]<6 or A.shape[1]<3:reject["small_subweb"]+=1;continue
        species=A.index.tolist();birds=A.columns.tolist()
        R=np.array([tab[x][0]=="radial" for x in species])
        if min(R.sum(),(~R).sum())<3:reject["insufficient_symmetry"]+=1;continue
        bb=rows.groupby("a").bill.agg(lambda s:float(s.iloc[0]) if s.nunique()==1 else np.nan)
        bill=bb.reindex(birds).to_numpy()
        if not np.isfinite(bill).all() or min(bill)<=0:reject["bill_conflict"]+=1;continue
        x=np.array([tab[s][1] for s in species])
        L=np.abs(x[:,None]-np.log2(bill)[None,:])
        A=A.to_numpy(dtype=bool);edges=np.argwhere(A)
        degree=(A.sum(0).copy(),A.sum(1).copy())
        n=len(edges);obs=metrics(A,L,R)
        moved=swap(edges,A,rng,20*n)
        if moved<max(10,n):reject["nonmixing_burnin"]+=1;continue
        sim=[];moves=[]
        for k in range(draws):
            moves.append(swap(edges,A,rng,5*n));sim.append(metrics(A,L,R))
        assert np.array_equal(A.sum(0),degree[0]) and np.array_equal(A.sum(1),degree[1])
        # Require successful switches over the retained chain; otherwise
        # uninformative margins cannot count as independent evidence.
        if np.mean(moves)<max(2,.1*n):reject["nonmixing_chain"]+=1;continue
        loc=rows[["lat","lon"]].drop_duplicates()
        if len(loc)!=1 or not np.isfinite(loc.to_numpy()).all():reject["location_unresolved"]+=1;continue
        lat,lon=loc.iloc[0]
        if not (-90<=lat<=90 and -180<=lon<=180):reject["invalid_location"]+=1;continue
        out.append({"block":(int(np.floor(lat/10)),int(np.floor(lon/10))),"region":"tropical" if abs(lat)<=23.5 else ("north" if lat>0 else "south"),"species":len(species),"birds":len(birds),"links":n,"radial":int(R.sum()),"obs":obs,"sim":np.array(sim),"moves":float(np.mean(moves))})
    return out,dict(reject)

def aggregate(rr):
    if not rr:return {"status":"no_eligible_networks"}
    keys=sorted({r["block"] for r in rr})
    a=np.stack([np.mean([r["obs"] for r in rr if r["block"]==b],axis=0) for b in keys])
    s=np.stack([np.mean([r["sim"] for r in rr if r["block"]==b],axis=0) for b in keys])
    o=a.mean(0);e=s.mean((0,1));res=o-e;draw=s.mean(0)
    rng=np.random.default_rng(SEED);n=len(keys)
    boot=np.stack([(a-s.mean(1))[rng.integers(n,size=n)].mean(0) for _ in range(999)])
    loo=np.stack([np.delete(a-s.mean(1),i,axis=0).mean(0) for i in range(n)]) if n>=2 else None
    pv=[min(1.,2*min((1+int(np.sum(draw[:,i]<=o[i]+1e-12)))/(len(draw)+1),(1+int(np.sum(draw[:,i]>=o[i]-1e-12)))/(len(draw)+1))) for i in range(3)]
    return {"networks":len(rr),"blocks":n,"gate":len(rr)>=20 and n>=10,"metrics":{k:{"observed":float(o[i]),"null_mean":float(e[i]),"residual":float(res[i]),"conditional_p":float(pv[i]),"block_residual_interval":np.quantile(boot[:,i],[.025,.975]).tolist() if n>=10 else None,"leave_block_range":[float(loo[:,i].min()),float(loo[:,i].max())] if loo is not None else None} for i,k in enumerate(("radial_mismatch","bilateral_mismatch","radial_minus_bilateral"))}}

def self_test():
    A=np.array([[1,0,1],[0,1,0],[1,1,0]],bool);E=np.argwhere(A);d=(A.sum(0).copy(),A.sum(1).copy())
    assert swap(E,A,np.random.default_rng(1),1000)>0
    assert np.array_equal(A.sum(0),d[0]) and np.array_equal(A.sum(1),d[1])
    assert np.array_equal(np.argwhere(A),np.array(sorted(map(tuple,E))))
    print("FD_HUMMINGBIRD_MATCHING_SELF_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args();self_test()
    if a.self_test:return
    d,source=fetch_network()
    d["p"]=d.Plant_species.map(binomial);d["a"]=d.Animal_species.map(binomial)
    d["y"]=pd.to_numeric(d.outcome,errors="coerce")
    d["bill"]=pd.to_numeric(d["Beak.Length_Culmen"],errors="coerce")
    d["lat"]=pd.to_numeric(d.Latitude,errors="coerce");d["lon"]=pd.to_numeric(d.Longitude,errors="coerce")
    assert d.y.isin([0,1]).all()
    audit={"rows":len(d),"networks":int(d.Network_name.nunique()),"plant_names":int(d.p.nunique()),"bird_names":int(d.a.nunique()),"positive_links":int(d.y.sum()),"unresolved_plant_rows":int(d.p.isna().sum()),"unresolved_bird_rows":int(d.a.isna().sum())}
    print("FD_HUMMINGBIRD_NETWORK_CENSUS",json.dumps(audit),flush=True)
    # Reject networks with unresolved names rather than treating them as absent.
    bad=set(d.loc[d.p.isna()|d.a.isna(),"Network_name"]);d=d[~d.Network_name.isin(bad)].copy()
    d.loc[~np.isfinite(d.bill)|(d.bill<=0),"bill"]=np.nan
    _,pins,nhm,wang=source_tables(include_values=True)
    sy=wang["unannotated_without_four_families"];ll=nhm["corolla_tube_range_cm"]
    results={}
    for scope in ("unannotated_without_four_families","all_publisher_binary"):
        for mode in ("midpoint","max"):
            sy=wang[scope]
            tab={s:(sy[s],float(np.log2(10)+(np.log2(ll[s][0])+np.log2(ll[s][1]))/2 if mode=="midpoint" else np.log2(10*ll[s][1]))) for s in sy.keys()&ll.keys()}
            rr,rej=panels(d,tab)
            key=scope+"__"+mode
            results[key]={"global_carriers":len(tab),"network_plant_overlap":len(set(d.p)&tab.keys()),"unresolved_networks":len(bad),"rejections":rej,"all":aggregate(rr),"tropical":aggregate([z for z in rr if z["region"]=="tropical"]),"north":aggregate([z for z in rr if z["region"]=="north"]),"south":aggregate([z for z in rr if z["region"]=="south"])}
            print("FD_HUMMINGBIRD_MATCHING_PANEL",key,json.dumps(results[key],allow_nan=False),flush=True)
    primary=results["unannotated_without_four_families__midpoint"]["all"];m=primary.get("metrics",{}).get("radial_minus_bilateral",{});c=m.get("block_residual_interval");l=m.get("leave_block_range")
    verdict=bool(primary.get("gate") and m.get("conditional_p",1)<.05 and c and l and ((c[0]>0 and l[0]>0) or (c[1]<0 and l[1]<0)))
    result={"source":source,"network_census":audit,"plant_source_pins":pins["source_pins"],"panels":results,"primary_supported":verdict,"warning":"Post-result independent-link calibration, not new island incidence or effective pollination. Observed versus unobserved binary links on documented local complete subwebs, filtered to exact measured NHM tube ranges (not effective nectar depth) and Wang symmetry. Culmen length not effective reach including tongue. Degree-preserving checkerboard Monte Carlo within each network, 399 draws, no formal chain mixing proof; geography block sensitivity, conditional on observed carrier subset and degree sequence. Primary all-network radial-minus-bilateral mismatch residual two-sided, >=20 networks and >=10 blocks, conditional p and consistent block interval/leave-block signs. Other scopes and region summaries descriptive; no island/mainland classification or claim of guild shift. Original source and repeated species overlap across networks."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False))
    print("FD_HUMMINGBIRD_MATCHING_COMPLETE",json.dumps({"panels":len(results),"primary_supported":verdict}),flush=True)
if __name__=="__main__":main()
