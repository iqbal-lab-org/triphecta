import filecmp
import logging
import os
import pytest
import subprocess
from unittest import mock

from triphecta import distances, tasks, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "tasks")


def test_pheno_constraints_template():
    options = mock.Mock()
    options.phenos_tsv = os.path.join(data_dir, "pheno_constraints_template.tsv")
    options.json_out = "tmp.tasks.pheno_constraints_template.json"
    utils.rm_rf(options.json_out)
    tasks.pheno_constraints_template.run(options)
    expect_file = os.path.join(data_dir, "pheno_constraints_template.json")
    assert filecmp.cmp(options.json_out, expect_file)
    os.unlink(options.json_out)


def test_pipeline(caplog):
    caplog.set_level(logging.INFO)
    vcf_names_file = "tmp.tasks.vcfs_to_names.tsv"
    mask_bed_file = os.path.join(data_dir, "mask.bed")
    distance_matrix_prefix = "tmp.distance_matrix"
    dist_matrix_file = f"{distance_matrix_prefix}.distance_matrix.tsv.gz"
    variant_counts_file = f"{distance_matrix_prefix}.variant_counts.tsv.gz"
    utils.rm_rf(vcf_names_file, dist_matrix_file, variant_counts_file)

    # ------------------ vcfs_to_names ----------------------------------------
    options = mock.Mock()
    options.file_of_vcf_filenames = os.path.join(data_dir, "vcfs.fofn")
    options.out_tsv = vcf_names_file
    options.threads = 1
    tasks.vcfs_to_names.run(options)
    expect = os.path.join(data_dir, "vcfs_to_names.expect.tsv")
    assert filecmp.cmp(options.out_tsv, expect, shallow=False)

    # ----------------- distance_matrix ---------------------------------------
    options = mock.Mock()
    options.method = "vcf"
    options.out = distance_matrix_prefix
    options.filenames_tsv = vcf_names_file
    options.threads = 1
    options.vcf_numeric_filter = None
    options.het_to_hom_key = None
    options.het_to_hom_cutoff = None
    options.mask_bed_file = mask_bed_file
    options.vcf_ignore_filter_pass = True
    expect_matrix_file = os.path.join(data_dir, "distance_matrix.tsv")
    expect_names, expect_distances = distances.load_distance_matrix_file(
        expect_matrix_file
    )
    tasks.distance_matrix.run(options)
    assert os.path.exists(dist_matrix_file)
    assert os.path.exists(variant_counts_file)
    got_names, got_distances = distances.load_distance_matrix_file(dist_matrix_file)
    assert got_names == expect_names
    assert got_distances == expect_distances

    # ----------------- tree --------------------------------------------------
    options = mock.Mock()
    options.distance_matrix = dist_matrix_file
    options.out = "tmp.tree.out"

    for method in "nj", "upgma":
        utils.rm_rf(options.out)
        options.method = method
        tasks.tree.run(options)
        assert os.path.exists(options.out)
        os.unlink(options.out)

    # ----------------- triples -----------------------------------------------
    options = mock.Mock()
    options.vcfs_tsv = vcf_names_file
    options.distance_matrix = dist_matrix_file
    options.var_counts_file = variant_counts_file
    options.phenos_tsv = os.path.join(data_dir, "phenos.tsv")
    options.pheno_constraints_json = os.path.join(data_dir, "pheno_constraint.json")
    options.out = "tmp.tasks.triples.out"
    utils.rm_rf("{options.out}.*")
    options.wanted_pheno = ["drug1,R"]
    options.top_n_genos = 5
    options.max_pheno_diffs = 1
    options.mask_bed_file = mask_bed_file
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
    os.unlink(dist_matrix_file)
    os.unlink(variant_counts_file)
    subprocess.check_output(f"rm -r {options.out}.triples", shell=True)
