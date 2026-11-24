#!/usr/bin/env python3
"""Audit provisional bee-bird candidate pairs and exact bee-ID butterfly bridge.
Descriptive diagnostics only: no certified bird geometry, guild completeness, or causal effect."""
import argparse,hashlib,io,json,math,os,pathlib,subprocess,zipfile
import numpy as np
import pandas as pd
PINS={
"bee_bird":(14901918629,"ed0cd97b3168d3e2618256ae9704aec409abbf5697c49d1dd1f81f72cb350a94","candidate_bee_bird_pairs.csv"),
"butterfly":(14896477401,"1d32946410e31f5df4999f488f98c307454835c82cdd1fe2c2a33720f2cd3472","butterfly_bee_positive_opportunity.csv")}
def fetch(kind):
    aid,digest,name=PINS[kind]
    url=f"https://api.github.com/repos/zuizui0223/fd/actions/artifacts/{aid}/zip"
    raw=subprocess.check_output(["curl","--fail","--silent","--show-error","--location","--max-time","90","--max-filesize","2000000","-H","Authorization: Bearer "+os.environ["GH_TOKEN"],"-H","Accept: application/vnd.github+json",url])
    assert hashlib.sha256(raw).hexdigest()==digest
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        paths=[x for x in z.namelist() if pathlib.PurePosixPath(x).name==name]
        assert len(paths)==1
        return pd.read_csv(z.open(paths[0]),dtype={"bee_id":str})
def slopes(d,dist,weight="group",resample=None):
    xs=[];ys=[];ng=0
    for g,v in d.groupby("group"):
        if len(v)<2 or v[dist].nunique()<2:continue
        x=np.column_stack([np.log2(1+v[dist].to_numpy(float)),np.log2(v.bee_area_km2.to_numpy(float))])
        y=v[["native_bees","nectar_ge50","nectar_ge20","nectar_gt0","non_nectar","all_birds"]].to_numpy(float)
        y=np.log1p(y)
        y=np.column_stack([y,y[:,1]-y[:,0],y[:,2]-y[:,0],y[:,3]-y[:,0],y[:,3]-y[:,4]])
        m=1 if resample is None else resample.get(g,0)
        if not m:continue
        w=math.sqrt(m/(len(v) if weight=="group" else 1))
        xs.append((x-x.mean(axis=0))*w);ys.append((y-y.mean(axis=0))*w);ng+=1
    if not xs:return None
    x=np.concatenate(xs);y=np.concatenate(ys)
    if np.linalg.matrix_rank(x)<2 or len(x)-ng<=2:return None
    names=["bees","nectar_ge50","nectar_ge20","nectar_gt0","non_nectar","all_birds","ge50_minus_bees","ge20_minus_bees","gt0_minus_bees","gt0_minus_non_nectar"]
    return dict(zip(names,np.linalg.lstsq(x,y,rcond=None)[0][0].tolist()))
def audit(d):
    out={};rng=np.random.default_rng(20261124)
    for label,v in [("all_candidates",d),("normalized_name_only",d[d.match_type=="normalized_name"]),("area_10pct",d[d.area_ratio.between(.9,1.1)]),("exclude_shared_boundary",d[~d.bird_name.isin(["St.Martin","Guadeloupe"])]),("strict_lat_bbox",d[d.bbox.map(lambda x:all(abs(float(t))<=23.44 for t in x.split(";")[1::2]))])]:
        v=v.copy()
        for dist in ("dist_km","dist_gt1","dist_gt1k"):
            v2=v[np.isfinite(v[dist])&(v[dist]>=0)]
            groups=sorted(g for g,w in v2.groupby("group") if len(w)>=2 and w[dist].nunique()>=2)
            result={"n":len(v2),"groups":int(v2.group.nunique()),"informative_groups":groups}
            for weight in ("group","island"):
                est=slopes(v2,dist,weight)
                drops={g:slopes(v2[v2.group!=g],dist,weight) for g in groups}
                result[weight]={"slopes":est,"leave_one_group":drops}
                if label=="all_candidates" and dist=="dist_km" and est:
                    vals=[]
                    for i in range(999):
                        a=rng.choice(groups,size=len(groups),replace=True)
                        q={g:int(np.sum(a==g)) for g in groups}
                        s=slopes(v2,dist,weight,q)
                        if s:vals.append(list(s.values()))
                    result[weight]["group_bootstrap_exploratory"]={k:np.quantile(np.array(vals)[:,j],[.025,.5,.975]).tolist() for j,k in enumerate(est)} if vals else None
            out[label+"_"+dist]=result
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    b=fetch("bee_bird");t=fetch("butterfly")
    assert len(b)==52 and b.bee_id.is_unique and len(t)==98 and t.bee_id.is_unique
    b["non_nectar"]=b.all_birds-b.nectar_gt0
    assert (b.non_nectar>=0).all()
    d=b[(b.group!="Ryukyu sensitivity")&(b.native_bees>0)].copy()
    assert len(d)==38 and d.group.nunique()==7
    joint=b.merge(t,on="bee_id",validate="one_to_one",suffixes=("_bird","_butterfly"))
    assert (joint.native_bees_bird==joint.native_bees_butterfly).all()
    bridge=[]
    for label,v in [("all",joint),("tropical_candidate",joint[joint.group!="Ryukyu sensitivity"])]:
        bridge.append({"scope":label,"candidate_units":len(v),"groups":sorted(v.group.unique()),"positive_butterfly":{k:int((v[k+"_records"]>0).sum()) for k in ("all","precision_le1km","no_gbifid_with_reference","both")},"all_three_positive":int(((v.all_records>0)&(v.nectar_gt0>0)&(v.native_bees_bird>0)).sum())})
    joint.to_csv(a.out/"three_guild_candidate_bridge.csv",index=False)
    r={"candidate_audit":audit(d),"bridge":bridge,"warning":"Provisional name/region/area match only. No independently certified bird geometry. Candidate groups are not sufficient for prespecified primary H1, and butterfly counts are positive records, not complete guild richness. Log1p count contrasts with non-nectar control are exploratory. Group resampling not a calibrated confidence interval."}
    (a.out/"summary.json").write_text(json.dumps(r,indent=2,allow_nan=False))
    print("THREE_GUILD_BRIDGE",json.dumps(bridge),flush=True)
    print("BIRD_BEE_CANDIDATE_SENSITIVITY",json.dumps(r["candidate_audit"],allow_nan=False),flush=True)
if __name__=="__main__":main()
