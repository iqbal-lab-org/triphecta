import copy
import numpy as np
import os

import pytest

from triphecta import vcf

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "vcf")


def test_vcf_line_to_variant_and_gt():
    fields = [
        "seq1",
        "42",
        "id",
        "C",
        "T,G",
        "42.42",
        "PASS",
        "INFO",
        "GT:FOO",
        "0/0:BAR",
    ]
    expect_variant = vcf.Variant(CHROM="seq1", POS=41, REF="C", ALTS=["T", "G"])
    got_gt, got_variant = vcf.vcf_line_to_variant_and_gt("\t".join(fields))
    assert got_gt == 0
    assert got_variant == expect_variant

    fields = [
        "seq2",
        "43",
        "id",
        "A",
        "C",
        "42.42",
        "FAIL",
        "INFO",
        "GT:FOO",
        "0/0:BAR",
    ]
    expect_variant = vcf.Variant(CHROM="seq2", POS=42, REF="A", ALTS=["C"])
    got_gt, got_variant = vcf.vcf_line_to_variant_and_gt("\t".join(fields))
    assert got_gt == "."
    assert got_variant == expect_variant

    fields[6] = "PASS"
    fields[8] = "X:Y"  # does not start with GT
    with pytest.raises(RuntimeError):
        vcf.vcf_line_to_variant_and_gt("\t".join(fields))

    fields[8] = "GT"
    fields[9] = "0/1"
    expect_variant = vcf.Variant(CHROM="seq2", POS=42, REF="A", ALTS=["C"])
    got_gt, got_variant = vcf.vcf_line_to_variant_and_gt("\t".join(fields))
    assert got_gt == "."
    assert got_variant == expect_variant

    fields[9] = "1/1"
    expect_variant = vcf.Variant(CHROM="seq2", POS=42, REF="A", ALTS=["C"])
    got_gt, got_variant = vcf.vcf_line_to_variant_and_gt("\t".join(fields))
    assert got_gt == 1
    assert got_variant == expect_variant


def test_load_variant_calls_from_vcf_file():
    infile = os.path.join(data_dir, "load_variants_from_vcf_file.vcf")
    expect_calls = [0, ".", 1, ".", "."]
    expect_variants = [
        vcf.Variant(CHROM="ref_42", POS=10, REF="C", ALTS=["G"]),
        vcf.Variant(CHROM="ref_42", POS=11, REF="A", ALTS=["C"]),
        vcf.Variant(CHROM="ref_43", POS=41, REF="T", ALTS=["A", "CT"]),
        vcf.Variant(CHROM="ref_43", POS=42, REF="T", ALTS=["C"]),
        vcf.Variant(CHROM="ref_43", POS=43, REF="T", ALTS=["G"]),
    ]
    got_calls, got_variants = vcf.load_variant_calls_from_vcf_file(infile)
    assert got_calls == expect_calls
    assert got_variants == expect_variants

    got_calls, got_variants = vcf.load_variant_calls_from_vcf_file(
        infile, expected_variants=expect_variants
    )
    assert got_calls == expect_calls
    assert got_variants == expect_variants

    expect_variants = copy.copy(got_variants)
    expect_variants[0] = vcf.Variant(CHROM="wrong_ref", POS=10, REF="C", ALTS=["G"])
    with pytest.raises(RuntimeError):
        vcf.load_variant_calls_from_vcf_file(infile, expected_variants=expect_variants)

    expect_variants = copy.copy(got_variants[:-1])
    with pytest.raises(RuntimeError):
        vcf.load_variant_calls_from_vcf_file(infile, expected_variants=expect_variants)


def test_load_vcf_file_for_distance_calc():
    infile = os.path.join(data_dir, "load_vcf_file_for_distance_calc.vcf")
    got_genos, got_counts = vcf.load_vcf_file_for_distance_calc(infile)
    expect_genos = np.array([0, 0, 2, 3, 4], dtype=np.uint16)
    np.testing.assert_array_equal(got_genos, expect_genos)
    expect_counts = vcf.VariantCounts(het=1, hom=3, null=1)
    assert got_counts == expect_counts

    got_genos, got_counts = vcf.load_vcf_file_for_distance_calc(
        infile, only_use_pass=False
    )
    expect_genos = np.array([1, 0, 2, 3, 4], dtype=np.uint16)
    np.testing.assert_array_equal(got_genos, expect_genos)
    expect_counts = vcf.VariantCounts(het=1, hom=4, null=0)
    assert got_counts == expect_counts

    got_genos, got_counts = vcf.load_vcf_file_for_distance_calc(
        infile, only_use_pass=False, numeric_filters={"GT_CONF": (True, 12)}
    )
    expect_genos = np.array([1, 0, 0, 0, 4], dtype=np.uint16)
    np.testing.assert_array_equal(got_genos, expect_genos)
    expect_counts = vcf.VariantCounts(het=1, hom=2, null=2)
    assert got_counts == expect_counts


def test_load_vcf_files_for_distance_calc():
    filenames = [
        os.path.join(data_dir, f"load_vcf_files_for_distance_calc.{i}.vcf")
        for i in (1, 2)
    ]

    def check_got_equal_expect(got, expect):
        assert len(got) == len(expect)
        for g, e in zip(got, expect):
            np.testing.assert_array_equal(g[0], e[0])
            assert g[1] == e[1]

    got = vcf.load_vcf_files_for_distance_calc(filenames, threads=2, only_use_pass=True)
    expect = [
        (np.array([0, 0, 2, 3, 4], dtype=np.uint16), vcf.VariantCounts(het=1, hom=3, null=1)),
        (np.array([0, 2, 2, 0, 2], dtype=np.uint16), vcf.VariantCounts(het=0, hom=3, null=2)),
    ]
    check_got_equal_expect(got, expect)

    got = vcf.load_vcf_files_for_distance_calc(
        filenames, threads=2, only_use_pass=False
    )
    expect = [
        (np.array([1, 0, 2, 3, 4], dtype=np.uint16), vcf.VariantCounts(het=1, hom=4, null=0)),
        (np.array([1, 2, 2, 3, 2], dtype=np.uint16), vcf.VariantCounts(het=0, hom=5, null=0)),
    ]
    check_got_equal_expect(got, expect)

    got = vcf.load_vcf_files_for_distance_calc(
        filenames,
        threads=2,
        only_use_pass=False,
        numeric_filters={"GT_CONF": (True, 12)},
    )
    expect = [
        (np.array([1, 0, 0, 0, 4], dtype=np.uint16), vcf.VariantCounts(het=1, hom=2, null=2)),
        (np.array([1, 2, 2, 0, 2], dtype=np.uint16), vcf.VariantCounts(het=0, hom=4, null=1)),
    ]
    check_got_equal_expect(got, expect)
