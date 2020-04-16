import os
import pytest
import subprocess

from triphecta import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_open_file():
    tmp_file = "tmp.open_file"
    subprocess.check_output(f"rm -f {tmp_file}", shell=True)
    with pytest.raises(OSError):
        with utils.open_file(tmp_file) as f:
            pass

    for tmp_file in "tmp.open_file", "tmp.open_file.gz":
        with utils.open_file(tmp_file, "w") as f:
            print("TEST", file=f)
            print("TEST2", file=f)
        assert os.path.exists(tmp_file)
        with utils.open_file(tmp_file) as f:
            lines = [x.rstrip() for x in f]
        assert lines == ["TEST", "TEST2"]
        os.unlink(tmp_file)


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
