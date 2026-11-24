#!/usr/bin/env python3
"""Actual LepTraits / butterfly incidence analysis. Trait filtering on RECORDS, not true colonization.
Fixed species-island incidence, observed trait masks and within-family/genus permutations.
"""
import argparse, collections, hashlib, io, json, math, os, pathlib, subprocess, zipfile
import numpy as np
import pandas as pd
from tropical_leptraits_intake import URL,BLOB,norm
from tropical_pollinator_source_preflight import get,blobsha

PIN=(14902243892,"99a72f6e6b5d72299917224a82d7d6e7cba17596cc47216ed5edc77c26dc3735")
MASKS=["all","precision_le1km","no_gbifid_with_reference","both"]
def load():
    aid,digest=PIN
    url=f"https://api.github.com/repos/zuizui0223/fd/actions/artifacts/{aid}/zip"
    raw=subprocess.check_output(["curl","--fail","--silent","--show-error","--location","--max-time","90","--max-filesize","2000000","-H","Authorization: Bearer "+os.environ["GH_TOKEN"],"-H","Accept: application/vnd.github+json",url])
    assert hashlib.sha256(raw).hexdigest()==digest
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        def read(name):
            p=[n for n in z.namelist() if pathlib.PurePosixPath(n).name==name];assert len(p)==1
            return pd.read_csv(z.open(p[0]),dtype={"bee_id":str})
        inc=read("butterfly_positive_incidence.csv"); geo=read("butterfly_bee_positive_opportunity.csv")
    assert len(geo)==98 and geo.bee_id.is_unique and len(inc)==18283
    inc["name"]=inc.final_species_name.map(norm)
    invalid=inc.name.isna().sum()
    inc=inc.dropna(subset=["name"]).drop_duplicates(["bee_id","source_mask","name"])
    normalized_count=int(inc.name.nunique())
    print("NORMALIZATION_REPLAY",json.dumps({"prior_logged_names":2313,"replayed_names":normalized_count,"invalid_rows":int(invalid),"deduplicated_rows":len(inc)}),flush=True)
    if normalized_count<1000: raise ValueError("unexpectedly few accepted binomials")
    raw=get(URL,maxbytes=5_000_000);assert blobsha(raw)==BLOB
    assert hashlib.sha256(raw).hexdigest()=="6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b"
    traits=pd.read_csv(io.BytesIO(raw),dtype=str)
    traits["name"]=traits.Species.map(norm)
    assert len(traits)==12448
    # Unambiguous duplicate consensus only. Never average discordant taxa or fill missing traits.
    cols=["name","Family","Genus","WS_L","WS_U","FW_L","FW_U","NumberOfHostplantFamilies","NumberOfHostplantAccounts"]
    tr=traits[cols].dropna(subset=["name"]).drop_duplicates()
    conflicting=set(tr.loc[tr.name.duplicated(False),"name"])
    tr=tr[~tr.name.isin(conflicting)].copy()
    assert tr.name.is_unique
    for col in cols[3:]:
        tr[col]=pd.to_numeric(tr[col],errors="coerce")
    def size(lo,hi):
        # Publication reports cm; require both source bounds and do not mix forewing/wingspan.
        a,b=tr[lo],tr[hi]; good=a.gt(0)&b.ge(a)&b.le(35)
        return np.log(((a+b)/2).where(good))
    tr["log_wingspan"]=size("WS_L","WS_U")
    tr["log_forewing"]=size("FW_L","FW_U")
    h=tr.NumberOfHostplantFamilies
    a=tr.NumberOfHostplantAccounts
    tr["log_host_families"]=np.log(h.where(h.gt(0)&h.le(200)&h.eq(np.floor(h))))
    tr["log_host_accounts"]=np.log1p(a.where(a.ge(1)))
    tr["family"]=tr.Family.fillna("").str.strip()
    tr["genus"]=tr.name.str.split().str[0]
    tr=tr[tr.family.ne("")]
    assert tr.name.is_unique
    coverage={"source_positive_rows":18283,"invalid_name_rows":int(invalid),"normalized_names":normalized_count,"prior_logged_names":2313,"exact_matches_before_duplicate_audit":int(inc.name.isin(traits.name.dropna()).groupby(inc.name).any().sum()),"prior_logged_exact_matches":1723,"conflicting_names":len(conflicting),"joined_unique_names":int(inc.name.isin(tr.name).groupby(inc.name).any().sum()),"trait_observed_distinct":{k:int(tr.loc[tr.name.isin(inc.name),k].notna().sum()) for k in ["log_wingspan","log_forewing","log_host_families","log_host_accounts"]},"units":"cm for wing measurements; log midpoint of positive lower/upper bounds. Host family count is documented breadth, not actual local host availability."}
    return inc,geo,tr,coverage
def group_residual(df,y,effort,weight="group"):
    xs=[];ys=[];ids=[]
    for g,v in df.groupby("archipelago"):
        if len(v)<3 or v.distance_km.nunique()<2:continue
        x=np.column_stack([np.log2(1+v.distance_km),np.log2(v.area_km2)]+([np.log2(v.records)] if effort else []))
        w=1/math.sqrt(len(v)) if weight=="group" else 1
        xs.append((x-x.mean(axis=0))*w);ys.append((y[v.index]-y[v.index].mean())*w);ids.extend(v.index)
    if not xs:return None
    x=np.concatenate(xs);yy=np.concatenate(ys)
    if np.linalg.matrix_rank(x)<x.shape[1] or len(yy)-len(xs)<=x.shape[1]:return None
    return float(np.linalg.lstsq(x,yy,rcond=None)[0][0])
def build_design(g,effort):
    xs=[];rowblocks=[];weight=[]
    for _,v in g.groupby("archipelago"):
        if len(v)<3 or v.distance_km.nunique()<2:continue
        x=np.column_stack([np.log2(1+v.distance_km),np.log2(v.area_km2)]+([np.log2(v.records)] if effort else []))
        w=1/math.sqrt(len(v))
        xs.append((x-x.mean(axis=0))*w);rowblocks.append(v.index.to_numpy());weight.append(w)
    if not xs:return None
    x=np.concatenate(xs)
    if np.linalg.matrix_rank(x)<x.shape[1] or len(x)-len(xs)<=x.shape[1]:return None
    coef=np.linalg.pinv(x)[0]
    # Map the coefficient to island means; centering each block.
    b=np.zeros(len(g));j=0
    for idx,w in zip(rowblocks,weight):
        v=coef[j:j+len(idx)]*w;b[idx]=v-v.mean();j+=len(idx)
    return b
def fit_endpoint(inc,geo,tr,mask,trait,minrec,mintrait,draws):
    ncol=mask+"_records"; g=geo[(geo[ncol]>=minrec)&(geo.area_km2>0)&(geo.distance_km>=0)].copy()
    g["records"]=g[ncol]; tab=inc[(inc.source_mask==mask)&inc.bee_id.isin(g.bee_id)].merge(tr[["name","family","genus",trait]],on="name",how="left",validate="many_to_one")
    cov=tab.groupby("bee_id")[trait].agg(["size","count"]).reset_index().rename(columns={"size":"recorded_names","count":"trait_names"})
    g=g.merge(cov,on="bee_id",how="inner",validate="one_to_one")
    g["coverage"]=g.trait_names/g.recorded_names
    # Require minimum observed traits and 50% coverage, preserve exact matched cohort for all permutations.
    g=g[(g.trait_names>=mintrait)&(g.coverage>=.5)].reset_index(drop=True)
    tab=tab[tab.bee_id.isin(g.bee_id)&tab[trait].notna()].copy()
    if not len(g):return {"n":0}
    species=tab[["name","family","genus",trait]].drop_duplicates().sort_values("name").reset_index(drop=True)
    assert species.name.is_unique
    row={v:i for i,v in enumerate(g.bee_id)}; col={v:i for i,v in enumerate(species.name)}
    matrix=np.zeros((len(g),len(species)))
    for r in tab.itertuples():matrix[row[r.bee_id],col[r.name]]=1
    matrix/=matrix.sum(axis=1,keepdims=True)
    x=species[trait].to_numpy(float);means=matrix@x
    out={"n":len(g),"groups":int(g.archipelago.nunique()),"informative_groups":int(sum(len(v)>=3 and v.distance_km.nunique()>=2 for _,v in g.groupby("archipelago"))),"species":len(species),"median_trait_coverage":float(g.coverage.median()),"results":{}}
    rng=np.random.default_rng(20261124)
    for effort in (False,True):
        b=build_design(g,effort)
        if b is None:continue
        coeff=b@matrix;obs=float(coeff@x)
        res={"observed_logtrait_slope_per_isolation_doubling":obs,"coverage_slope":float(b@g.coverage.to_numpy(float)),"permutation":{}}
        for tax in ("family","genus"):
            groups=[v.index.to_numpy() for _,v in species.groupby(tax) if len(v)>1 and v[trait].nunique()>1]
            if not groups:continue
            # Fixed incidence and observation mask; distinct species permuted coherently across all islands.
            sims=np.empty(draws)
            for i in range(draws):
                xp=x.copy()
                for ids in groups:xp[ids]=rng.permutation(x[ids])
                sims[i]=coeff@xp
            m=float(np.median(sims));p=(1+int((np.abs(sims-m)>=abs(obs-m)).sum()))/(draws+1)
            res["permutation"][tax]={"movable_species":sum(map(len,groups)),"blocks":len(groups),"null_q025_q50_q975":np.quantile(sims,[.025,.5,.975]).tolist(),"observed_minus_null_median":obs-m,"two_sided_centered_p":p}
        drops={}
        for group in sorted(g.archipelago.unique()):
            maskg=g.archipelago!=group;sub=g[maskg].reset_index(drop=True)
            z=build_design(sub,effort);drops[group]=float(z@means[maskg.to_numpy()]) if z is not None else None
        res["leave_one_archipelago"]=drops;out["results"]["effort_adjusted" if effort else "raw"]=res
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--draws",type=int,default=999);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    inc,g,tr,coverage=load();out={}
    for mask in MASKS:
        for trait in ("log_wingspan","log_forewing","log_host_families","log_host_accounts"):
            for minrec in (10,50):
                out[mask+"_"+trait+"_"+str(minrec)]=fit_endpoint(inc,g,tr,mask,trait,minrec,5,a.draws)
    result={"coverage":coverage,"panels":out,"warning":"Exploratory many related source record tests; no corrected significance or causal dispersal. Broad Asian/Pacific box is not strict tropical; record-count adjustment is diagnostic only. Trait values constant across islands; family/genus shuffles retain incidence and trait documentation masks but do not prove phylogenetic independence. Host breadth and wing size are not local resources or measured flight."}
    (a.out/"butterfly_trait_filtering.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    print("BUTTERFLY_TRAIT_COVERAGE",json.dumps(coverage),flush=True)
    print("BUTTERFLY_TRAIT_FILTERING",json.dumps(out),flush=True)
if __name__=="__main__":main()
