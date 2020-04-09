import filecmp
import logging
import os
import pytest
import subprocess
from unittest import mock

from triphecta import tasks

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "tasks")


def test_pheno_constraints_template():
    options = mock.Mock()
    options.phenos_tsv = os.path.join(data_dir, "pheno_constraints_template.tsv")
    options.json_out = "tmp.tasks.pheno_constraints_template.json"
    subprocess.check_output(f"rm -f {options.json_out}", shell=True)
    tasks.pheno_constraints_template.run(options)
    expect_file = os.path.join(data_dir, "pheno_constraints_template.json")
    assert filecmp.cmp(options.json_out, expect_file)
    os.unlink(options.json_out)


def test_pipeline(caplog):
    caplog.set_level(logging.INFO)
    vcf_names_file = "tmp.tasks.vcfs_to_names.tsv"
    mask_bed_file = os.path.join(data_dir, "mask.bed")
    distances_file = "tmp.distances.pickle"
    subprocess.check_output(f"rm -rf {vcf_names_file} {distances_file}", shell=True)

    # ------------------ vcfs_to_names ----------------------------------------
    options = mock.Mock()
    options.file_of_vcf_filenames = os.path.join(data_dir, "vcfs.fofn")
    options.out_tsv = vcf_names_file
    options.threads = 1
    tasks.vcfs_to_names.run(options)
    expect = os.path.join(data_dir, "vcfs_to_names.expect.tsv")
    assert filecmp.cmp(options.out_tsv, expect, shallow=False)

    # ------------------ distances --------------------------------------------
    options = mock.Mock()
    options.method = "vcf"
    options.outfile = distances_file
    options.filenames_tsv = vcf_names_file
    options.threads = 1
    options.vcf_numeric_filter = None
    options.het_to_hom_key = None
    options.het_to_hom_cutoff = None
    options.mask_bed_file = mask_bed_file
    options.vcf_ignore_filter_pass = True
    tasks.distances.run(options)
    assert os.path.exists(distances_file)

    # ----------------- triples -----------------------------------------------
    options = mock.Mock()
    options.vcfs_tsv = vcf_names_file
    options.distances_file = distances_file
    options.phenos_tsv = os.path.join(data_dir, "phenos.tsv")
    options.pheno_constraints_json = os.path.join(data_dir, "pheno_constraint.json")
    options.out = "tmp.tasks.triples.out"
    options.wanted_pheno = ["drug1,R"]
    options.top_n_genos = 5
    options.max_pheno_diffs = 1
    options.mask_file = mask_bed_file
    tasks.triples.run(options)

    got_triple_ids_tsv = f"{options.out}.triple_ids.tsv"
    expect_triple_ids_tsv = os.path.join(data_dir, "triples.triple_ids.tsv")
    assert filecmp.cmp(got_triple_ids_tsv, expect_triple_ids_tsv, shallow=False)
    os.unlink(got_triple_ids_tsv)

    got_variants_tsv = f"{options.out}.variants.tsv"
    expect_variants_tsv = os.path.join(data_dir, "triples.variants.tsv")
    assert filecmp.cmp(got_variants_tsv, expect_variants_tsv, shallow=False)
    os.unlink(got_variants_tsv)

    for i in (1, 2, 3):
        expect_tsv = os.path.join(data_dir, "triples", f"{i}.tsv")
        got_tsv = os.path.join(f"{options.out}.triples", f"{i}.tsv")
        assert filecmp.cmp(got_tsv, expect_tsv, shallow=False)

    os.unlink(vcf_names_file)
    os.unlink(distances_file)
    subprocess.check_output(f"rm -r {options.out}.triples", shell=True)
