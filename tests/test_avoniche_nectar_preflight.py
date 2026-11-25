"""Synthetic AVONICHE file checks; no GIFT IDs or bird polygons used."""
import pandas as pd
import pytest
from fd.avoniche_nectar_preflight import inspect_public_birdlife_csv


def source(rows):
    return pd.DataFrame(rows).to_csv(index=False).encode("utf-8")


def rows():
    return [dict(Scientific_Name="Bird first", Pa=0, Pt=0, Ne=20, Se=0,
                 Fr=30, In=50, Ap=0, Vt=0, Ca=0),
            dict(Scientific_Name="Bird second", Pa=0, Pt=0, Ne=0, Se=0,
                 Fr=50, In=50, Ap=0, Vt=0, Ca=0)]


def test_percentage_scale_is_certified_by_all_nine_diets():
    receipt, weights = inspect_public_birdlife_csv(source(rows()), expected_rows=2)
    assert receipt["diet_scale"] == 100
    assert receipt["n_recorded_nectar_use_positive"] == 1
    assert weights["nectar_fraction"].tolist() == [.2, 0.]
    assert receipt["range_polygon_join_completed"] is False


def test_fraction_scale_is_certified_by_nine_diets():
    data = [{k: (v / 100 if k in {"Pa","Pt","Ne","Se","Fr","In","Ap","Vt","Ca"} else v)
             for k, v in r.items()} for r in rows()]
    receipt, _ = inspect_public_birdlife_csv(source(data), expected_rows=2)
    assert receipt["diet_scale"] == 1


def test_fails_when_nectar_ratio_unverifiable():
    d = rows()
    for r in d: del r["Fr"]
    with pytest.raises(ValueError, match="without independent nine-diet"):
        inspect_public_birdlife_csv(source(d), expected_rows=2)


def test_fails_on_duplicate_species_and_missing_diet():
    d = rows()
    d[1]["Scientific_Name"] = d[0]["Scientific_Name"]
    with pytest.raises(ValueError, match="duplicated"):
        inspect_public_birdlife_csv(source(d), expected_rows=2)
    d = rows()
    d[0]["Ne"] = None
    with pytest.raises(ValueError, match="missing"):
        inspect_public_birdlife_csv(source(d), expected_rows=2)


def test_refuses_wrong_source_population_or_ambiguous_columns():
    with pytest.raises(ValueError, match="expected 10999"):
        inspect_public_birdlife_csv(source(rows()))
    d = rows()
    for r in d: r["Nectar_raw"] = r.pop("Ne")
    with pytest.raises(ValueError, match="Unverified nectar diet schema"):
        inspect_public_birdlife_csv(source(d), expected_rows=2)
