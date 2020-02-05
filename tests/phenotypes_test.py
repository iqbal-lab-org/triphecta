import os

import pytest

from triphecta import phenotypes

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "phenotypes")


def test_convert_one_variable_string():
    assert not phenotypes.Phenotypes._convert_one_variable_string("f")
    assert not phenotypes.Phenotypes._convert_one_variable_string("s")
    assert not phenotypes.Phenotypes._convert_one_variable_string("susceptible")
    assert phenotypes.Phenotypes._convert_one_variable_string("r")
    assert phenotypes.Phenotypes._convert_one_variable_string("resistant")
    assert phenotypes.Phenotypes._convert_one_variable_string("true")
    assert phenotypes.Phenotypes._convert_one_variable_string("") is None
    assert phenotypes.Phenotypes._convert_one_variable_string(" ") is None
    assert phenotypes.Phenotypes._convert_one_variable_string(".") is None
    assert phenotypes.Phenotypes._convert_one_variable_string("None") is None
    assert phenotypes.Phenotypes._convert_one_variable_string("null") is None
    with pytest.raises(ValueError):
        phenotypes.Phenotypes._convert_one_variable_string("what is this")
    with pytest.raises(ValueError):
        phenotypes.Phenotypes._convert_one_variable_string("1.2.3")


def test_load_phenotypes_tsv_file():
    infile = os.path.join(data_dir, "load_phenotype_file.tsv")
    got_pheno, got_types = phenotypes.Phenotypes._load_phenotypes_tsv_file(infile)
    expected_pheno = {
        "s1": {"pheno1": 1, "pheno2": True},
        "s2": {"pheno1": 2, "pheno2": False},
    }
    assert got_pheno == expected_pheno
    expected_types = {"pheno1": {float}, "pheno2": {bool}}
    assert got_types == expected_types


def test_get_pheno_types():
    types = {"p1": {float}, "p2": {bool}, "p3": {float, type(None)}}
    got_all, got_bools = phenotypes.Phenotypes._get_pheno_types(types)
    expect_all = {"p1": float, "p2": bool, "p3": float}
    expect_bools = {"p2"}
    assert got_all == expect_all
    assert got_bools == expect_bools
    types["p4"] = {float, bool}
    with pytest.raises(RuntimeError):
        phenotypes.Phenotypes._get_pheno_types(types)
    types["p4"] = {type(None)}
    expect_all["p4"] = type(None)
    got_all, got_bools = phenotypes.Phenotypes._get_pheno_types(types)
    assert got_all == expect_all
    assert got_bools == expect_bools
