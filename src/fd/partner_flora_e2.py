"""Source-agnostic fixed-incidence partner / plant E1-E2 diagnostic.
Frozen extension contract docs/PARTNER_E2_EXTENSION_CONTRACT_20261125.md.
No range or checklist absence inference; caller certifies data/source eligibility.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import sparse
from .nectarivore_flora_test import add_controls, group_index, plant_block_expectation, within_fit, permute_within_blocks, CONTROLS
from .range_map_nectarivore import diet_counts

def _design(inc, labels, d, mode):
    c = diet_counts(inc, labels)
    if mode=="diet":
        x=np.column_stack([np.log1p(c[:,2]),np.log1p(c[:,1]),np.log1p(c[:,0]),d[list(CONTROLS)].to_numpy(float)])
        contrast=np.array([1.,-1.,0.,0.,0.,0.,0.])
    elif mode=="fruit":
        x=np.column_stack([np.log1p(c[:,1]+c[:,2]),np.log1p(c[:,0]),d[list(CONTROLS)].to_numpy(float)])
        contrast=np.array([1.,0.,0.,0.,0.,0.])
    elif mode=="bombus":
        if (c[:,1]>0).any():raise ValueError("Bombus labels must be 0 or 2")
        x=np.column_stack([(c[:,2]>0).astype(float),np.log1p(c[:,0]),d[list(CONTROLS)].to_numpy(float)])
        contrast=np.array([1.,0.,0.,0.,0.,0.])
    else:raise ValueError("invalid partner mode")
    return x,contrast

def _block_ids(tab,cols):
    return pd.factorize(pd.MultiIndex.from_frame(tab[list(cols)]),sort=True)[0]

def _movable(labels,blocks):
    return sum(int(np.sum(blocks==g)) for g in np.unique(blocks) if len(np.unique(labels[blocks==g]))>1)

def _p(draw,obs):
    return (1+int(np.sum(np.asarray(draw)>=obs-1e-12)))/(len(draw)+1)

def holm_fixed(pvalues, ids=("NECTAR_BIRD","BOMB_NORTH","NECTAR_BAT")):
    """Unsupported slots = 1; reject nonfamily keys and invalid p values."""
    if not set(pvalues)<=set(ids):raise ValueError("unknown primary family")
    p=np.array([pvalues.get(k,1.) if pvalues.get(k) is not None else 1. for k in ids],float)
    if not np.isfinite(p).all() or (p<0).any() or (p>1).any():raise ValueError("invalid p")
    out=np.empty(len(p));running=0.
    for j,i in enumerate(np.argsort(p)):
        running=max(running,(len(p)-j)*p[i]);out[i]=min(1.,running)
    return dict(zip(ids,map(float,out)))

def run_partner_test(islands, plants, incidence, animal_labels, animal_genus=None, mode="diet",
                     *, min_groups=30, min_islands=200, n_boot=1999, n_perm=1999, seed=20261125):
    """islands: one row per eligible island, island/group/area_km2/dist_km/lat.
    plants: complete fixed regional pool of scored island-species rows (island/species/genus/grade/trait).
    incidence rows align to islands; only qualified source-supported units may be supplied.
    Counts, E2 and permutations never use an unobserved animal as an inferred absence.
    """
    if min_groups<1 or min_islands<1 or n_boot<1 or n_perm<1:raise ValueError("invalid gate")
    d=islands.copy().reset_index(drop=True)
    if d.island.isna().any() or d.island.duplicated().any() or d.group.isna().any() or d.group.astype(str).str.strip().eq("").any():raise ValueError("invalid island/group key")
    inc=sparse.csr_matrix(incidence)
    labels=np.asarray(animal_labels)
    if inc.shape!=(len(d),len(labels)) or not np.isin(labels,[0,1,2]).all() or not np.isin(inc.data,[0,1]).all():raise ValueError("bad incidence or animal labels")
    if plants[["species","genus","grade","trait"]].isna().any().any():raise ValueError("unknown plant state")
    shares=plant_block_expectation(plants)
    d=d.join(shares,on="island",validate="one_to_one")
    if d.n.isna().any() or (d.n<5).any() or (d.area_km2<1).any() or (d.dist_km<0).any() or (d.lat.abs()>90).any():raise ValueError("unqualified island")
    d=add_controls(d)
    # Keep the fixed regional plant pool, but permute animals in the actually present eligible pool only.
    present=np.asarray(inc.sum(axis=0)).ravel()>0
    inc=inc[:,present];labels=labels[present]
    ag=None if animal_genus is None else np.asarray(animal_genus)[present]
    if ag is not None and pd.isna(ag).any():raise ValueError("unknown animal taxonomic block")
    X,c=_design(inc,labels,d,mode)
    if not np.isfinite(X).all():raise ValueError("nonfinite design")
    groups=group_index(d.group)
    informative=sum(len(ix)>=2 and np.ptp(X[ix,0])>0 for ix in groups)
    ninfo=sum(len(ix) for ix in groups if len(ix)>=2)
    out={"mode":mode,"n_islands":len(d),"groups":len(groups),"informative_groups":informative,"within_group_islands":ninfo,
         "gate_groups":min_groups,"gate_islands":min_islands,"animal_species_present":int(present.sum()),"status":"insufficient_support",
         "p_primary":None}
    if informative<min_groups or len(d)<min_islands:return out
    dof=sum(len(ix)-1 for ix in groups)-X.shape[1]
    if dof<4:return out|{"status":"insufficient_rank"}
    tab=plants.drop_duplicates("species").reset_index(drop=True)
    ids={s:i for i,s in enumerate(tab.species)}
    idx={s:i for i,s in enumerate(d.island)}
    pl=plants[plants.island.isin(idx)]
    ri=pl.island.map(idx).to_numpy(int);ci=pl.species.map(ids).to_numpy(int)
    A=sparse.csr_matrix((np.ones(len(pl)),(ri,ci)),shape=(len(d),len(tab)))
    A=sparse.diags(1/d.n.to_numpy(float))@A
    V=tab.trait.to_numpy(float);blocks=_block_ids(tab,("genus","grade"))
    means=pd.Series(V).groupby(blocks).transform("mean").to_numpy()
    expected=np.asarray(A@means).ravel();y=np.asarray(A@V).ravel()
    if not np.allclose(y,d.share) or not np.allclose(expected,d.expected):raise AssertionError("expectation mismatch")
    resid=y-expected
    def est(x,v,mult=None):
        b=within_fit(x,v,groups,mult)
        return None if b is None else float(c@b)
    if sum(np.var(resid[ix]) for ix in groups)<1e-12:return out|{"status":"not_estimable"}
    e1,e2=est(X,y),est(X,resid)
    if e1 is None or e2 is None:return out|{"status":"insufficient_rank"}
    out.update(E1=e1,E2=e2,E_expected=e1-e2,dof=dof,plant_species=len(tab),plant_movable_species=_movable(V,blocks),status="estimated")
    rng=np.random.default_rng(seed);draws=[]
    for _ in range(n_boot):
        mult=np.bincount(rng.integers(0,len(groups),len(groups)),minlength=len(groups))
        v=est(X,resid,mult)
        if v is not None:draws.append(v)
    out["bootstrap_valid"]=len(draws);out["E2_interval"]=np.quantile(draws,[.025,.975]).tolist() if draws else None
    if len(draws)<.9*n_boot or out["plant_movable_species"]==0:return out|{"status":"insufficient_null_support"}
    if mode!="bombus" and (ag is None or _movable(labels,ag)==0):return out|{"status":"insufficient_animal_null_support","animal_movable_species":0}
    out["animal_movable_species"]=None if mode=="bombus" else _movable(labels,ag)
    pdraw=[];adraw=[];jdraw=[]
    for _ in range(n_perm):
        shuffled=permute_within_blocks(V,blocks,rng)
        ry=np.asarray(A@shuffled).ravel()-expected
        v=est(X,ry)
        if v is not None:pdraw.append(v)
        if mode!="bombus":
            lab=permute_within_blocks(labels,ag,rng)
            nx,_=_design(inc,lab,d,mode)
            av,jv=est(nx,resid),est(nx,ry)
            if av is not None and jv is not None:adraw.append(av);jdraw.append(jv)
    out["permutation_valid"]={"plant":len(pdraw),"animal":len(adraw),"joint":len(jdraw)}
    if len(pdraw)<.9*n_perm or (mode!="bombus" and len(adraw)<.9*n_perm):return out|{"status":"insufficient_null_support"}
    ps={"plant":_p(pdraw,e2)}
    if mode!="bombus":ps.update(animal_genus=_p(adraw,e2),joint=_p(jdraw,e2))
    out["permutation_p"]=ps
    out["p_primary"]=max(ps.values()) if out["E2_interval"][0]>0 else 1.
    out["status"]="ok"
    return out
