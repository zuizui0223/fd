#!/usr/bin/env python3
"""Official-source, post-acquisition breeding-system calibration.

A measured AF index is reproductive performance under pollinator exclusion,
not a realized selfing rate and not proof of sexual autonomous selfing.
Compare only exact normalized binomials and same-row, same-metric SC/AF.
Data_source labels and GIFT reference ID sets are provenance strata, not
verified independent primary studies. No island or own-project inputs.
"""
import argparse,collections,hashlib,io,json,pathlib,re,time,urllib.parse
import numpy as np
import pandas as pd
from tropical_pollinator_source_preflight import get
from external_breeding_system_source_audit import CLOUD,cloud_media_url,sha

BOOK="PMC5095580.1/ncomms13313-s2.xlsx"
BOOK_SHA="cb17009ae65c4f6cb2d28b21239714782600d0904effa77c5635249d28e62fc0"
META_SHA="09db64b758571279494350bc558c9b1c51515aee6239cf3273e2ee058a5ce410"
API="https://gift.uni-goettingen.de/api/extended/index3.1.php"
PINS={"guild":("3.6.2","66cfa0398c3222131e642f083e567651efc1719e04b15a18c1462911988e6b75",5852),
      "self":("3.1.1","59d44b73c6267fb6bb08a677b58395c5480392415a50b5eccbfd9e111407fd17",3530)}
METRICS={"FS":("Self-compatibility_index_FS","Autofertility_index_FS"),
         "SFL":("Self-compatibility_index_SFL","Autofertility_index_SFL")}
GROUPS=(("wind","absent"),("wind","present"),("insect","absent"),("insect","present"))
SEED=20261125
DRAWS=1999

def name(s):
    # Formatting normalization only. No author stripping, synonym
    # inference, infraspecies collapse or fuzzy match.
    if pd.isna(s):return None
    x=" ".join(str(s).strip().replace("_"," ").split()).casefold()
    return x if re.fullmatch(r"[a-z][a-z-]+ [a-z][a-z-]+",x) else None

def clean(s):
    return "" if pd.isna(s) else str(s).strip()

def fetch(query,**params):
    url=API+"?"+urllib.parse.urlencode({"query":query,**params})
    b=get(url,maxbytes=80000000)
    rows=json.loads(b)
    assert isinstance(rows,list) and all(isinstance(r,dict) for r in rows)
    return rows,{"url":url,"sha256":sha(b),"rows":len(rows)}

def refs(v):
    return tuple(sorted({x for x in re.split(r"[,;|\s]+",str(v or "")) if x.isdigit()}))

def traits():
    result={};pins={}
    for k,(tid,pin,n) in PINS.items():
        rows,p=fetch("traits",traitid=tid,biasref=0,biasderiv=0,startat=0,limit=10000)
        assert p["sha256"]==pin
        d={}
        for r in rows:
            i=str(r["work_ID"]);assert i not in d
            v=str(r["trait_value"]).strip().casefold()
            a=r.get("agreement")
            if a not in (None,"") and float(a)<.66:continue
            if v in ("","none","na","nan","null"):continue
            d[i]=(v,refs(r.get("references")))
        assert len(d)==n,(k,len(d),n)
        result[k]=d;pins[k]=p
    return result,pins

def book():
    prefix="PMC5095580.1/"
    raw=get(CLOUD+prefix+prefix.rstrip("/")+".json",maxbytes=2000000)
    assert sha(raw)==META_SHA
    m=json.loads(raw)
    assert m["doi"]=="10.1038/ncomms13313" and m["pmcid"]=="PMC5095580"
    assert m["is_pmc_openaccess"] is True and m["is_retracted"] is False
    assert str(m["license_code"]).lower()=="cc by"
    matches=[u for u in m["media_urls"] if urllib.parse.urlsplit(u).path.endswith("/ncomms13313-s2.xlsx")]
    assert len(matches)==1
    url,k,md5=cloud_media_url(matches[0],prefix);assert k==BOOK
    b=get(url,maxbytes=1000000)
    assert hashlib.md5(b).hexdigest()==md5 and sha(b)==BOOK_SHA
    xl=pd.ExcelFile(io.BytesIO(b),engine="openpyxl")
    d=xl.parse("Razanajatovo_etal_data")
    remarks=xl.parse("Remarks")
    assert len(d)==1890 and {"Species","Family","Data_source"}<=set(d.columns)
    assert all(c in d.columns for pair in METRICS.values() for c in pair)
    return d,{"workbook_url":url,"sha256":BOOK_SHA,"publisher_md5":md5,"metadata_sha256":META_SHA,
              "rows":len(d),"remarks":remarks.fillna("").astype(str).to_dict(orient="records")}

def names_index(wanted):
    # Official GIFT R package retrieves the species endpoint in six
    # 100,000-row pages. Track all work_ID matches for each wanted name.
    d=collections.defaultdict(set);pins=[]
    for offset in range(0,600000,100000):
        rows,p=fetch("species",startat=offset)
        assert len(rows)<=100000
        pins.append(p)
        for r in rows:
            x=name(r.get("work_species"))
            if x in wanted:d[x].add(str(r["work_ID"]))
        if len(rows)<100000:break
        time.sleep(.2)
    else:
        raise ValueError("species pagination cap reached; incomplete ambiguity audit")
    return d,pins

def cohort(d,metric,index,T):
    sc,af=METRICS[metric]
    x=d.copy()
    x["name"]=x.Species.map(name)
    x["sc"]=pd.to_numeric(x[sc],errors="coerce")
    x["af"]=pd.to_numeric(x[af],errors="coerce")
    ok=x.sc.between(0,1)&x.af.between(0,1)
    counts={"all_rows":len(x),"sc_in_unit_interval":int(x.sc.between(0,1).sum()),
            "af_in_unit_interval":int(x.af.between(0,1).sum()),"same_row_pair_valid":int(ok.sum()),
            "same_row_pair_valid_named":int((ok&x.name.notna()).sum()),
            "nonmissing_sc_invalid":int((x[sc].notna()&~x.sc.between(0,1)).sum()),
            "nonmissing_af_invalid":int((x[af].notna()&~x.af.between(0,1)).sum())}
    rows=[];reasons=collections.Counter()
    for n,g in x[ok&x.name.notna()].groupby("name",sort=True):
        ids=index.get(n,set())
        if len(ids)!=1:reasons["no_unique_gift_work_id"]+=1;continue
        i=next(iter(ids))
        if i not in T["guild"] or i not in T["self"]:reasons["both_gift_traits_missing"]+=1;continue
        gv,gr=T["guild"][i];sv,sr=T["self"][i]
        if (gv,sv) not in GROUPS:reasons["outside_frozen_insect_wind_self_cells"]+=1;continue
        if not gr or not sr:reasons["gift_reference_ids_missing"]+=1;continue
        fam=sorted({clean(v).casefold() for v in g.Family if clean(v)})
        if len(fam)!=1:reasons["missing_or_conflicting_family"]+=1;continue
        src=tuple(sorted({clean(v) for v in g.Data_source if clean(v)}))
        if not src:reasons["workbook_source_label_missing"]+=1;continue
        rows.append({"name":n,"work_id":i,"genus":n.split()[0],"family":fam[0],
            "sc":float(g.sc.mean()),"af":float(g.af.mean()),"paired_rows":len(g),
            "guild":gv,"self":sv,"source_pattern":json.dumps(src),
            "gift_pattern":json.dumps([gr,sr])})
    counts["distinct_paired_normalized_names"]=int(x.loc[ok,"name"].nunique())
    counts["attrition"]=dict(reasons)
    counts["eligible_species"]=len(rows)
    assert len({r["name"] for r in rows})==len({r["work_id"] for r in rows})==len(rows)
    return pd.DataFrame(rows),counts

def describe(d):
    if d.empty:return {"species":0,"cells":[]}
    cells=[]
    for g,s in GROUPS:
        a=d[(d.guild==g)&(d.self==s)]
        cells.append({"guild":g,"gift_self":s,"species":len(a),"genera":a.genus.nunique(),
            "families":a.family.nunique(),"SC_mean":float(a.sc.mean()) if len(a) else None,
            "AF_mean":float(a.af.mean()) if len(a) else None,
            "SC_min_max":[float(a.sc.min()),float(a.sc.max())] if len(a) else None,
            "AF_zero":int((a.af==0).sum()),"AF_one":int((a.af==1).sum())})
    return {"species":len(d),"genera":d.genus.nunique(),"families":d.family.nunique(),
        "workbook_source_patterns":d.source_pattern.nunique(),"gift_reference_patterns":d.gift_pattern.nunique(),
        "workbook_source_pattern_counts":{k:int(v) for k,v in d.source_pattern.value_counts().head(12).items()},
        "gift_reference_pattern_counts":{k:int(v) for k,v in d.gift_pattern.value_counts().head(12).items()},
        "duplicate_paired_rows":int(d.paired_rows.sum()-len(d)),"cells":cells}

def fit(d,kind):
    if d.empty:return {"status":"not_estimable","support":describe(d)}
    # Common SC support across four cells; chosen from predictor values
    # alone. Avoid extrapolating one guild into an unobserved SC range.
    bounds=[]
    for g,s in GROUPS:
        z=d[(d.guild==g)&(d.self==s)].sc
        if len(z)<20:return {"status":"insufficient_cell_support","support":describe(d)}
        bounds.append((float(z.min()),float(z.max())))
    lo=max(x[0] for x in bounds);hi=min(x[1] for x in bounds)
    if hi-lo<.1:return {"status":"insufficient_common_SC_range","SC_common_range":[lo,hi],"support":describe(d)}
    a=d[d.sc.between(lo,hi)].copy()
    center=(lo+hi)/2
    a["stratum"]="all" if kind=="raw" else a.family if kind=="family" else a.family+"|"+a.source_pattern+"|"+a.gift_pattern
    # Remove strata that cannot identify any conditional contrast. The
    # strict provenance+family variant requires both guilds and both
    # coded self states in each retained stratum.
    if kind!="raw":
        keep=[]
        for k,g in a.groupby("stratum"):
            cell=g.groupby(["guild","self"]).size()
            if all(cell.get(z,0)>=2 for z in GROUPS):keep.append(k)
        a=a[a.stratum.isin(keep)].copy()
    support=describe(a);support["strata"]=int(a.stratum.nunique()) if len(a) else 0
    if len(a)==0 or min(c["species"] for c in support["cells"])<20 or min(c["genera"] for c in support["cells"])<10 or support["genera"]<30 or (kind!="raw" and support["strata"]<5):
        return {"status":"insufficient_matched_support","support":support,"SC_common_range":[lo,hi]}
    s=(a.self=="present").to_numpy(float);g=(a.guild=="insect").to_numpy(float);sc=a.sc.to_numpy(float)-center
    X=np.column_stack([sc,s,g,s*g]);y=a.af.to_numpy(float)
    group=pd.factorize(a.stratum,sort=True)[0];gen=pd.factorize(a.genus,sort=True)[0]
    def solve(w):
        counts=np.bincount(group,weights=w)
        z=np.column_stack([X,y])
        means=np.column_stack([np.bincount(group,weights=w*z[:,j])/counts for j in range(z.shape[1])])
        zz=(z-means[group])*np.sqrt(w[:,None])
        beta,_,rank,sv=np.linalg.lstsq(zz[:,:4],zz[:,4],rcond=None)
        if rank!=4 or sv[-1]<=1e-10:return None
        return np.array([beta[1],beta[1]+beta[3],beta[3],beta[2],beta[2]+beta[3],beta[0]])
    point=solve(np.ones(len(a)))
    if point is None:return {"status":"rank_deficient","support":support}
    rng=np.random.default_rng(SEED);draw=[]
    for i in range(DRAWS):
        w=np.maximum(rng.exponential(size=int(gen.max())+1),1e-12)[gen]
        v=solve(w)
        if v is not None:draw.append(v)
    return {"status":"exploratory_estimable","support":support,"SC_common_range":[lo,hi],
        "SC_center":center,"metric_order":["coded_self_AF_association_wind","coded_self_AF_association_insect",
        "coded_self_by_guild_interaction","insect_minus_wind_when_coded_absent","insect_minus_wind_when_coded_present","SC_slope"],
        "point":point.tolist(),"valid_draws":len(draw),
        "genus_weight_sensitivity_interval":np.quantile(draw,[.025,.975],axis=0).T.tolist() if len(draw)>=.9*DRAWS else None,
        "warning":"Linear conditional descriptive contrast, shared positive genus weights, not calibrated sampling CI. Workbooks do not identify independent experiments per row. No causal guild or reproductive mechanism."}

def self_test():
    assert name("Abelia_chinensis")=="abelia chinensis"
    assert name("Abelia chinensis R.Br.") is None and name("Abelia chinensis subsp. x") is None
    # Same-row requirement: do not pair SC and AF on different rows.
    d=pd.DataFrame({"Species":["Aaa_aaa","Aaa_aaa","Bbb_bbb"],"Family":["F"]*3,"Data_source":["1"]*3,
         METRICS["FS"][0]:[.5,None,.8],METRICS["FS"][1]:[None,.7,.4]})
    idx={"aaa aaa":{"1"},"bbb bbb":{"2"}}
    T={"guild":{"1":("insect",("1",)),"2":("wind",("1",))},"self":{"1":("present",("2",)),"2":("absent",("2",))}}
    a,c=cohort(d,"FS",idx,T);assert len(a)==1 and a.iloc[0]["sc"]==.8
    assert fit(a,"raw")["status"]=="insufficient_cell_support"
    print("EXTERNAL_MEASURED_CALIBRATION_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    d,source=book()
    wanted={n for s in d.Species if (n:=name(s)) is not None}
    names_audit={"raw_rows":len(d),"nonmissing_name_rows":int(d.Species.notna().sum()),
        "normalized_binomial_rows":int(d.Species.map(name).notna().sum()),"distinct_normalized_binomials":len(wanted),
        "underscore_name_rows":int(d.Species.astype(str).str.contains("_",regex=False).sum()),
        "unparsed_examples":sorted({clean(x) for x in d.Species if clean(x) and name(x) is None})[:12]}
    print("EXTERNAL_MEASURED_NAMES",json.dumps(names_audit,ensure_ascii=False),flush=True)
    T,trait_pins=traits();index,species_pins=names_index(wanted)
    names_audit["unique_gift_name_matches"]=sum(len(v)==1 for v in index.values())
    names_audit["ambiguous_gift_name_matches"]=sum(len(v)>1 for v in index.values())
    print("EXTERNAL_MEASURED_NAME_MATCH",json.dumps(names_audit,ensure_ascii=False),flush=True)
    out={"source":source,"name_audit":names_audit,"trait_pins":trait_pins,"species_pins":species_pins,"panels":{},
         "warning":"Source-calibration only. AF indices from pollinator exclusion relative to hand pollination can include apomixis. SC and AF on same workbook row and metric, averaged across duplicate species rows. Exact normalized names only; no synonym inference. GIFT coded self fertilization not necessarily AF, coded syndromes not necessarily observed pollination. Source labels and reference ID sets are not established primary study independence. No island effect or causal evolutionary inference."}
    for metric in METRICS:
        x,census=cohort(d,metric,index,T)
        z={"census":census,"all_matched":describe(x),"fits":{k:fit(x,k) for k in ("raw","family","family_and_source_patterns")}}
        out["panels"][metric]=z
        print("EXTERNAL_MEASURED_CALIBRATION_PANEL",metric,json.dumps(z,ensure_ascii=False,allow_nan=False),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False))
    print("EXTERNAL_MEASURED_CALIBRATION_COMPLETE",json.dumps({"panels":len(out["panels"]),"source_sha256":BOOK_SHA}),flush=True)
if __name__=="__main__":main()
