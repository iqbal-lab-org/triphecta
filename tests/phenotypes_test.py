import filecmp
import os
import subprocess

import pytest

from triphecta import phenotype_compare, phenotypes

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "phenotypes")


def test_convert_one_variable_string():
    assert not phenotypes.Phenotypes.convert_one_variable_string("f")
    assert not phenotypes.Phenotypes.convert_one_variable_string("s")
    assert not phenotypes.Phenotypes.convert_one_variable_string("susceptible")
    assert phenotypes.Phenotypes.convert_one_variable_string("r")
    assert phenotypes.Phenotypes.convert_one_variable_string("resistant")
    assert phenotypes.Phenotypes.convert_one_variable_string("true")
    assert phenotypes.Phenotypes.convert_one_variable_string("") is None
    assert phenotypes.Phenotypes.convert_one_variable_string(" ") is None
    assert phenotypes.Phenotypes.convert_one_variable_string(".") is None
    assert phenotypes.Phenotypes.convert_one_variable_string("None") is None
    assert phenotypes.Phenotypes.convert_one_variable_string("null") is None
    with pytest.raises(ValueError):
        phenotypes.Phenotypes.convert_one_variable_string("what is this")
    with pytest.raises(ValueError):
        phenotypes.Phenotypes.convert_one_variable_string("1.2.3")


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


def test_write_template_constraints_json():
    phenos = phenotypes.Phenotypes(
        os.path.join(data_dir, "write_template_constraints_json.tsv")
    )
    tmp_json = "tmp.phenos.write_template_constraints.json"
    subprocess.check_output(f"rm -f {tmp_json}", shell=True)
    phenos.write_template_constraints_json(tmp_json)
    expect_json = os.path.join(data_dir, "write_template_constraints_json.json")
    assert filecmp.cmp(tmp_json, expect_json, shallow=True)
    os.unlink(tmp_json)


def test_find_matching_cases():
    phenos = phenotypes.Phenotypes(os.path.join(data_dir, "find_matching_cases.tsv"))
    constraints = {"pheno2": {"method": "equal", "must_be_same": True, "params": {},}}
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    got = phenos.find_matching_cases({"pheno2": "R"}, pheno_compare)
    assert got == ["s1"]
    got = phenos.find_matching_cases({"pheno2": "T"}, pheno_compare)
    assert got == ["s1"]
    got = phenos.find_matching_cases({"pheno2": "S"}, pheno_compare)
    assert got == ["s2", "s3"]

    constraints = {"pheno3": {"method": "equal", "must_be_same": False, "params": {},}}
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    got = phenos.find_matching_cases({"pheno3": "R"}, pheno_compare)
    assert got == ["s1", "s2"]
    got = phenos.find_matching_cases({"pheno3": "S"}, pheno_compare)
    assert got == []

    constraints = {
        "pheno1": {
            "method": "abs_distance",
            "must_be_same": False,
            "params": {"max_dist": 5},
        },
        "pheno2": {"method": "equal", "must_be_same": False, "params": {},},
    }
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    got = phenos.find_matching_cases({"pheno1": 1}, pheno_compare)
    assert got == ["s1", "s2"]
    got = phenos.find_matching_cases({"pheno1": 200}, pheno_compare)
    assert got == ["s3"]
    got = phenos.find_matching_cases({"pheno1": 100}, pheno_compare)
    assert got == []
    got = phenos.find_matching_cases({"pheno1": 100, "pheno2": "R"}, pheno_compare)
    assert got == []
    got = phenos.find_matching_cases({"pheno1": 1, "pheno2": "R"}, pheno_compare)
    assert got == ["s1"]
    got = phenos.find_matching_cases({"pheno1": 1, "pheno2": "S"}, pheno_compare)
    assert got == ["s2"]
    got = phenos.find_matching_cases({"pheno1": 199, "pheno2": "S"}, pheno_compare)
    assert got == ["s3"]
