"""Public DoPI v2 flower-visitor source audit, not a pollination-syndrome labeler.

The source (Balfour et al. 2022 and DoPI v2 2026, DOI
10.5281/zenodo.22876365) contains European POSITIVE plant-visitor edges.
It does not establish true absences, pollen-deposition effectiveness or
island-specific guild transitions. GloBI indexes DoPI, so de-duplicate
study provenance against GloBI before considering independent validation.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from hashlib import sha256
import csv
import json
from pathlib import Path
from urllib.request import Request,urlopen

RECORD_ID=22876365
DOI="10.5281/zenodo.22876365"
API=f"https://zenodo.org/api/records/{RECORD_ID}"
BEE_FAMILIES={"apidae","andrenidae","halictidae","colletidae",
              "megachilidae","melittidae","stenotritidae"}
BUTTERFLY_FAMILIES={"nymphalidae","papilionidae","pieridae","lycaenidae",
                    "hesperiidae","riodinidae"}
MOTH_FAMILIES={"sphingidae","noctuidae","erebidae","geometridae",
               "crambidae","pyralidae","saturniidae","bombycidae",
               "notodontidae","lasiocampidae"}
BIRD_ORDERS={"apodiformes","passeriformes","coliiformes"}


def identify_visitor_group(order:str|None,family:str|None,species:str|None)->str:
    o=(order or "").strip().lower()
    f=(family or "").strip().lower()
    s=(species or "").strip().lower()
    if s.startswith("bombus ") and f in BEE_FAMILIES:return "Bombus"
    if f in BEE_FAMILIES:return "other_bee"
    if f=="sphingidae":return "Sphingidae_hawkmoth"
    if f in BUTTERFLY_FAMILIES:return "butterfly_family"
    if f in MOTH_FAMILIES:return "other_identified_moth_family"
    if f=="syrphidae":return "Syrphidae_hoverfly"
    if o=="diptera":return "other_Diptera"
    if o=="lepidoptera":return "unresolved_Lepidoptera_family"
    if o in BIRD_ORDERS:return "bird_order_not_effective_pollinator"
    if o=="hymenoptera":return "other_Hymenoptera"
    if o=="coleoptera":return "Coleoptera"
    return "other_or_unresolved_visitor"


def source_metadata(record:dict)->dict:
    # Fail closed: DOI is the version-specific public release, not a mutable
    # current website or a different GloBI-ingested 2022 DoPI source.
    md=record.get("metadata",{})
    title=str(md.get("title",""))
    identifiers={str(x).strip().lower() for x in
                 (md.get("doi"),record.get("doi"),md.get("prereserve_doi",{}).get("doi")
                  if isinstance(md.get("prereserve_doi"),dict) else None)}
    if "dopi" not in title.lower() or DOI not in identifiers:
        raise ValueError("metadata is not immutable DoPI v2 release")
    files=record.get("files",[])
    matches=[f for f in files if
             str(f.get("key","")).lower() in {
                 "dopi_dataset.csv","dopi_dataset.csv.gz","dopi_dataset.zip"}]
    if len(matches)!=1:
        raise ValueError("expected single original DoPI_dataset original public release file")
    f=matches[0]
    url=f.get("links",{}).get("self")
    if not isinstance(url,str) or not url.startswith("https://"):
        raise ValueError("no authentic Zenodo DoPI file download URL")
    return dict(release_id=RECORD_ID,doi=DOI,
                title=title,version=md.get("version"),
                license=md.get("license"),
                dataset_file=f["key"],size_bytes=f.get("size"),
                checksum=f.get("checksum"),
                file_url=url)


def _norm(value):
    return str(value or "").strip()


def audit_positive_visitor_rows(rows, *,fields:set[str])->dict:
    minimum={"plant_species","pollinator_species",
             "pollinator_family","pollinator_order"}
    if not minimum<=fields:
        raise ValueError("DoPI taxonomic species, family and order fields absent")
    pubs=("publication" in fields or "doi" in fields)
    if not pubs:
        raise ValueError("DoPI independent original-study identity absent")
    by_guild=Counter()
    by_pub=defaultdict(set)
    by_plant=defaultdict(set)
    by_plant_pub_guild=defaultdict(set)
    counts=Counter()
    sample_effort=Counter()
    methods=Counter()
    years=set()
    for row in rows:
        counts["raw_rows"]+=1
        plant=_norm(row.get("plant_species"))
        animal=_norm(row.get("pollinator_species"))
        source=_norm(row.get("doi")) or _norm(row.get("publication"))
        if not plant or not animal or not source:
            counts["unusable_plant_animal_or_study"]+=1
            continue
        guild=identify_visitor_group(row.get("pollinator_order"),
                                     row.get("pollinator_family"),animal)
        counts["positive_taxon_study_rows"]+=1
        by_guild[guild]+=1
        by_pub[guild].add(source)
        by_plant[guild].add(plant)
        by_plant_pub_guild[guild].add((source,plant))
        sample_effort[(_norm(row.get("pollinator_survey")).lower(),
                       _norm(row.get("plant_survey")).lower())]+=1
        methods[_norm(row.get("methodology")).lower()]+=1
        y=_norm(row.get("year"))
        if y.isdecimal() and 1800<=int(y)<=2100:
            years.add(int(y))
    return dict(
        n_original_raw_records=int(counts["raw_rows"]),
        n_positive_taxon_study_rows=int(counts["positive_taxon_study_rows"]),
        n_rows_without_both_species_and_study=int(counts["unusable_plant_animal_or_study"]),
        n_studies_by_positive_guild={g:len(by_pub[g]) for g in sorted(by_pub)},
        n_plants_by_positive_guild={g:len(by_plant[g]) for g in sorted(by_plant)},
        n_study_plant_positive_units={g:len(by_plant_pub_guild[g]) for g in sorted(by_plant_pub_guild)},
        n_rows_by_positive_guild=dict(sorted(by_guild.items())),
        n_original_survey_scope_codes=len(sample_effort),
        top_survey_scope_labels=[[list(k),v] for k,v in sample_effort.most_common(12)],
        top_methods=[[k,v] for k,v in methods.most_common(12)],
        year_min=min(years) if years else None,
        year_max=max(years) if years else None,
        flags=dict(no_zero_absence_derived=True,
                   distinct_studies_not_proven_independent_of_GloBI=True,
                   visit_records_not_effective_pollination=True,
                   Europe_scope_not_global_tropical_bird_calibration=True,
                   do_not_use_study_plant_guild_positive_as_negative_labels=True))


def _get_json(url:str)->dict:
    request=Request(url,headers={"User-Agent":"fd-botanical-audit/0.1 (research; citation DOI:10.5281/zenodo.22876365)"})
    with urlopen(request,timeout=50) as resp:
        return json.load(resp)


def _download_to(path:Path,url:str,max_bytes:int=1500*1024*1024)->str:
    path.parent.mkdir(parents=True,exist_ok=True)
    sha=sha256()
    n=0
    request=Request(url,headers={"User-Agent":"fd-botanical-audit/0.1"})
    with urlopen(request,timeout=90) as resp,path.open("wb") as f:
        while True:
            block=resp.read(1<<20)
            if not block:break
            n+=len(block)
            if n>max_bytes:raise ValueError("Zenodo CSV unexpectedly large")
            f.write(block)
            sha.update(block)
    if n<100:raise ValueError("empty DoPI original source download")
    return sha.hexdigest()


def main(output:Path,workdir:Path):
    import gzip
    import zipfile
    metadata=_get_json(API)
    packet=source_metadata(metadata)
    rawfile=workdir/packet["dataset_file"]
    disk_sha=_download_to(rawfile,packet["file_url"])
    key=rawfile.name.lower()
    if key.endswith(".zip"):
        with zipfile.ZipFile(rawfile) as z:
            matches=[k for k in z.namelist() if Path(k).name.lower()=="dopi_dataset.csv"]
            if len(matches)!=1:
                raise ValueError("source ZIP missing original DoPI_dataset.csv")
            with z.open(matches[0]) as f:
                import io
                with io.TextIOWrapper(f,encoding="utf-8-sig",errors="replace",newline="") as handle:
                    reader=csv.DictReader(handle)
                    audit=audit_positive_visitor_rows(reader,fields=set(reader.fieldnames or []))
    else:
        opener=gzip.open if key.endswith(".gz") else open
        with opener(rawfile,"rt",encoding="utf-8-sig",errors="replace",newline="") as handle:
            reader=csv.DictReader(handle)
            audit=audit_positive_visitor_rows(reader,fields=set(reader.fieldnames or []))
    out=dict(status="PUBLIC_DOPI_V2_POSITIVE_VISITOR_PRECALIBRATION_AUDIT",
             original_file_sha256=disk_sha,
             original_metadata={k:v for k,v in packet.items() if k!="file_url"},
             **audit)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print(json.dumps(dict(status=out["status"],
                          n_original_raw_records=out["n_original_raw_records"],
                          n_positive_guild_groups=len(out["n_rows_by_positive_guild"]),
                          n_distinct_study_guilds=sum(out["n_studies_by_positive_guild"].values()))))


if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--workdir",type=Path,required=True)
    a=ap.parse_args()
    main(a.output,a.workdir)
