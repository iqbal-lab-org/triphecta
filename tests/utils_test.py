import os

import pytest

from triphecta import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_load_file_of_vcf_filenames():
    infile = os.path.join(data_dir, "load_file_of_vcf_filenames.tsv")
    vcf_files = ["sample1.vcf", "sample2.vcf"]
    for f in vcf_files:
        if os.path.exists(f):
            os.unlink(f)
    expect = {"sample1": "sample1.vcf", "sample2": "sample2.vcf"}
    with pytest.raises(FileNotFoundError):
        got = utils.load_file_of_vcf_filenames(infile)

    for f in vcf_files:
        with open(f, "w") as f:
            pass
    got = utils.load_file_of_vcf_filenames(infile)
    for f in vcf_files:
        os.unlink(f)
    assert got == expect


def test_command_line_filter_list_to_dict():
    filters_list = ["A:min:10", "B:max:20"]
    expect = {"A": (True, 10.0), "B": (False, 20.0)}
    got = utils.command_line_filter_list_to_dict(filters_list)
    assert got == expect


def test_command_line_wanted_phenos_to_dict():
    pheno_list = ["Drug1,r", "Drug2,42"]
    got = utils.command_line_wanted_phenos_to_dict(pheno_list)
    expect = {"Drug1": True, "Drug2": 42.0}
    assert got == expect
