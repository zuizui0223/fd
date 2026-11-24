#!/usr/bin/env python3
"""Sensitivity analysis of original source-derived public artifacts.
Observed record richness is NOT true guild richness. Bird checklist occupancy is NOT colonization.
"""
import argparse,hashlib,io,json,math,pathlib,os,subprocess,zipfile
from collections import defaultdict
import numpy as np
import pandas as pd

PINS={
"bird":(14896670354,"f9a948a35c7c89f1863953052cb1b4b94036811dbfe09ddeefd001c7efdbd53a","bird_hwi_occupancy.csv"),
"butterfly":(14896477401,"1d32946410e31f5df4999f488f98c307454835c82cdd1fe2c2a33720f2cd3472","butterfly_bee_positive_opportunity.csv")}
def load_artifact(kind):
    aid,expected,filename=PINS[kind]
    url=f"https://api.github.com/repos/zuizui0223/fd/actions/artifacts/{aid}/zip"
    raw=subprocess.check_output(["curl","--fail","--silent","--show-error","--location","--max-time","90","--max-filesize","2000000","-H","Authorization: Bearer "+os.environ["GH_TOKEN"],"-H","Accept: application/vnd.github+json",url])
    assert hashlib.sha256(raw).hexdigest()==expected,"artifact digest mismatch"
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        p=[n for n in z.namelist() if pathlib.PurePosixPath(n).name==filename]
        assert len(p)==1
        d=pd.read_csv(z.open(p[0]))
    return d

def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x,(np.bool_,)):return bool(x)
    return x
def quantile(v):
    v=[x for x in v if x is not None and np.isfinite(x)]
    return {"draws":len(v),"q025_q50_q975":np.quantile(v,[.025,.5,.975]).tolist() if v else None}

AQUATIC={"anatidae","podicipedidae","gaviidae","spheniscidae","diomedeidae","procellariidae","hydrobatidae","oceanitidae","phaethontidae","fregatidae","sulidae","phalacrocoracidae","anhingidae","pelecanidae","ardeidae","threskiornithidae","ciconiidae","phoenicopteridae","laridae","stercorariidae","alcidae","charadriidae","scolopacidae","haematopodidae","recurvirostridae","burhinidae","glareolidae","jacanidae","dromadidae","rallidae","heliornithidae","chionidae","ibidorhynchidae","rostratulidae","pluvianellidae"}
def cells_for(d,split=False,min_n=3):
    out=[]
    for (group,family),v in d.groupby(["group","family"]):
        if split:
            a=v.loc[v.nectar];b=v.loc[~v.nectar]
            if min(len(a),len(b))<min_n:continue
            if min(a.hwi.nunique(),b.hwi.nunique(),a.occupancy.nunique(),b.occupancy.nunique())<2:continue
            out.append((a.index.to_numpy(),b.index.to_numpy()))
        else:
            if len(v)<min_n or min(v.hwi.nunique(),v.occupancy.nunique())<2:continue
            out.append(v.index.to_numpy())
    return out
def moments(x,y,w,idx):
    w=w[idx];keep=w>0
    if keep.sum()<2:return None
    w=w[keep];xx=x[idx][keep];yy=y[idx][keep]
    w=w/w.sum();xx=xx-np.dot(w,xx);yy=yy-np.dot(w,yy)
    den=np.dot(w,xx*xx);num=np.dot(w,xx*yy)
    if den<1e-10:return None
    return num,den

def slope(x,y,w,cells,groups,group_mult=None):
    by=defaultdict(list)
    for idx in cells:
        m=moments(x,y,w,idx)
        if m is not None:by[groups[idx[0]]].append(m)
    if not by:return None
    # each active group contributes its mean cell moments; no unequal-cell group dominance
    a=np.asarray([np.mean(v,axis=0) for v in by.values()])
    m=np.array([1 if group_mult is None else group_mult.get(g,0) for g in by])
    num,den=(a*m[:,None]).sum(axis=0)
    return float(num/den*10) if den>0 else None
def split_slope(x,y,w,cells,groups,group_mult=None):
    am=[];bm=[]
    for a,b in cells:
        if moments(x,y,w,a) is not None and moments(x,y,w,b) is not None:
            am.append(a);bm.append(b)
    if not am:return None
    a=slope(x,y,w,am,groups,group_mult);b=slope(x,y,w,bm,groups,group_mult)
    if a is None or b is None:return None
    return a-b
def bird_audit(d,draws=1999):
    assert len(d)==820 and not d.duplicated(["group","species"]).any()
    d=d[d.group!="Ryukyu sensitivity"].copy()
    assert len(d)==721 and d.species.nunique()==497
    d["nectar"]=d.diet_nect>0
    results={}
    rng=np.random.default_rng(20261124)
    for scope,df in [("all",d),("non_aquatic_family",d[~d.family.str.lower().isin(AQUATIC)])]:
        for label,sel,split in [("nectar_gt0",df[df.nectar],False),("nectar_ge20",df[df.diet_nect>=20],False),("non_nectar",df[~df.nectar],False),("paired_family_interaction",df,True)]:
            sel=sel.reset_index(drop=True);n=len(sel)
            c=cells_for(sel,split=split,min_n=2 if split else 3)
            groups=sel.group.to_numpy();x=sel.hwi.to_numpy(float);y=sel.occupancy.to_numpy(float)
            gi=sorted(set(groups));sp=sorted(set(sel.species));gm={g:i for i,g in enumerate(gi)};sm={s:i for i,s in enumerate(sp)}
            gidx=np.array([gm[g] for g in groups]);sidx=np.array([sm[s] for s in sel.species]);w=np.ones(n)
            fx=split_slope if split else slope
            support_idx=np.unique(np.concatenate([np.concatenate(v) for v in c] if split else c)) if c else np.array([],dtype=int)
            support_groups=sorted(set(groups[support_idx]));support_species=len(set(sel.species.iloc[support_idx]))
            observed=fx(x,y,w,c,groups)
            if observed is None:
                results[scope+"_"+label]={"rows":n,"support_cells":len(c),"supported_groups":support_groups,"support_species":support_species,"slope":None}
                continue
            drops={g:fx(x,y,(groups!=g).astype(float),c,groups) for g in support_groups}
            # Independent resampling of archipelago and distinct species identities. A small number of groups makes this exploratory, not calibrated CI.
            vals=[];sg=[]
            for _ in range(draws):
                a=np.bincount(rng.integers(0,len(gi),len(gi)),minlength=len(gi))
                b=np.bincount(rng.integers(0,len(sp),len(sp)),minlength=len(sp))
                vals.append(fx(x,y,a[gidx]*b[sidx],c,groups,dict(zip(gi,a))));sg.append(fx(x,y,b[sidx],c,groups))
            # Conditional distinct-species HWI shuffle within family and genus separately.
            perm={}
            if not split:
                ss=sel.drop_duplicates("species")
                orig=dict(zip(ss.species,ss.hwi));species=np.array(sp);xsp=np.array([orig[z] for z in sp])
                assert all(sel.groupby("species").hwi.nunique()==1)
                for level in ("family","genus"):
                    tax=ss.family if level=="family" else ss.species.str.split().str[0]
                    taxmap=dict(zip(ss.species,tax))
                    blocks=defaultdict(list)
                    for i,z in enumerate(sp):blocks[taxmap[z]].append(i)
                    blocks=[np.array(v) for v in blocks.values() if len(v)>1 and np.ptp(xsp[v])>0]
                    v=[]
                    for _ in range(draws):
                        xp=xsp.copy()
                        for ix in blocks:xp[ix]=rng.permutation(xsp[ix])
                        v.append(fx(xp[sidx],y,w,c,groups))
                    valid=[a for a in v if a is not None]
                    perm[level]={"movable_species":sum(len(b) for b in blocks),"blocks":len(blocks),"null":quantile(valid),"one_sided_positive_p":(1+sum(a>=observed for a in valid))/(1+len(valid)) if valid else None}
            results[scope+"_"+label]={"rows":n,"support_cells":len(c),"supported_groups":support_groups,"support_species":support_species,"slope":observed,"leave_one_group":drops,"two_way_group_species_resampling":quantile(vals),"species_only_resampling":quantile(sg),"permutation":perm}
    return results

def within_group_regression(d,adjust_effort=False,equal_group=True):
    if len(d)<5:return None
    # paired log recorded counts; intercept removed by demeaning within original named archipelago.
    xs=[];ys=[];ng=0
    for _,v in d.groupby("archipelago"):
        if len(v)<2 or v.distance_km.nunique()<2:continue
        a=np.log2(1+v.distance_km.to_numpy(float))
        b=np.log2(v.area_km2.to_numpy(float))
        cols=[a,b]
        if adjust_effort:cols.append(np.log2(v.records.to_numpy(float)))
        x=np.column_stack(cols);y=np.log(v.butterfly.to_numpy(float))-np.log(v.native_bees.to_numpy(float))
        wt=math.sqrt(len(v)) if equal_group else 1
        xs.append((x-x.mean(0))/wt);ys.append((y-y.mean())/wt);ng+=1
    if not xs:return None
    x=np.concatenate(xs);y=np.concatenate(ys)
    if np.linalg.matrix_rank(x)<x.shape[1] or len(y)-ng<=x.shape[1]:return None
    return float(np.linalg.lstsq(x,y,rcond=None)[0][0])
def butterfly_audit(d):
    assert len(d)==98 and d.bee_id.nunique()==98
    endpoints=("all","precision_le1km","no_gbifid_with_reference","both")
    expected=[88,74,48,1]
    assert [int((d[k+"_records"]>0).sum()) for k in endpoints]==expected
    results={}
    # Full original bounding-box frame (includes extra-tropical / outside-source islands), do not call this strictly tropical.
    for endpoint in endpoints:
        for m in (1,10,50):
            v=d.loc[(d[endpoint+"_records"]>=m)&(d[endpoint+"_recorded_species"]>0)&(d.native_bees>0)&(d.area_km2>0)&(d.distance_km>=0)].copy()
            v["records"]=v[endpoint+"_records"];v["butterfly"]=v[endpoint+"_recorded_species"]
            support={str(g):len(q) for g,q in v.groupby("archipelago") if len(q)>=2 and q.distance_km.nunique()>=2}
            groups=sorted(support)
            res={"n_positive":len(v),"all_groups":int(v.archipelago.nunique()),"within_group_support":support}
            for effort in (False,True):
                k="effort_adjusted" if effort else "unadjusted"
                res[k]={"equal_group_slope":within_group_regression(v,effort),"island_equal_slope":within_group_regression(v,effort,False),"leave_one_group":{g:within_group_regression(v[v.archipelago!=g],effort) for g in groups}}
            results[endpoint+"_"+str(m)]=res
    # EXACT common positive support across all and high precision and non-GBIF-id/reference subsets (not precision AND non-GBIF).
    for m in (10,50):
        v=d.loc[np.logical_and.reduce([d[k+"_records"].ge(m) for k in endpoints[:3]]) & (d.native_bees>0)&(d.area_km2>0)&(d.distance_km>=0)].copy()
        res={"n":len(v),"groups":int(v.archipelago.nunique()),"endpoints":{}}
        for k in endpoints[:3]:
            q=v.assign(records=v[k+"_records"],butterfly=v[k+"_recorded_species"])
            res["endpoints"][k]={"raw":within_group_regression(q),"effort_adjusted":within_group_regression(q,True)}
        results["common_three_"+str(m)]=res
    return results

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--draws",type=int,default=1999);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    bird=bird_audit(load_artifact("bird"),a.draws);butterfly=butterfly_audit(load_artifact("butterfly"))
    out=clean({"pins":PINS,"bird":bird,"butterfly_bee":butterfly,"warning":"Exploratory source record/realized occupancy diagnostics only; family/genus permutations are not phylogenetic causal tests. No certified differential dispersal, pollination or GIFT floral linkage. Bird group/species resampling range not calibrated inferential interval. Butterfly effort adjustment is a sensitivity, not identification."})
    (a.out/"sensitivity.json").write_text(json.dumps(out,indent=2,allow_nan=False))
    print("BIRD_SENSITIVITY",json.dumps(out["bird"],allow_nan=False),flush=True)
    print("BUTTERFLY_BEE_SENSITIVITY",json.dumps(out["butterfly_bee"],allow_nan=False),flush=True)
if __name__=="__main__":main()
