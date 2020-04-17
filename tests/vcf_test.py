import copy
import filecmp
import numpy as np
import os

import pytest

from triphecta import variant_counts, vcf

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
    assert got_gt == {0}
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
    assert got_gt == None
    assert got_variant == expect_variant

    fields[6] = "PASS"
    fields[8] = "X:Y"  # does not start with GT
    with pytest.raises(RuntimeError):
        vcf.vcf_line_to_variant_and_gt("\t".join(fields))

    fields[8] = "GT"
    fields[9] = "0/1"
    expect_variant = vcf.Variant(CHROM="seq2", POS=42, REF="A", ALTS=["C"])
    got_gt, got_variant = vcf.vcf_line_to_variant_and_gt("\t".join(fields))
    assert got_gt == {0, 1}
    assert got_variant == expect_variant

    fields[9] = "1/1"
    expect_variant = vcf.Variant(CHROM="seq2", POS=42, REF="A", ALTS=["C"])
    got_gt, got_variant = vcf.vcf_line_to_variant_and_gt("\t".join(fields))
    assert got_gt == {1}
    assert got_variant == expect_variant


def test_load_variant_calls_from_vcf_file():
    infile = os.path.join(data_dir, "load_variants_from_vcf_file.vcf")
    expect_calls = [{0}, {0, 1}, {1}, None, None]
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


def test_convert_het_to_hom():
    genos = {"0", "1"}
    info = {"COV": "9,90,1"}
    assert vcf._convert_het_to_hom(genos, info, "wrong_key", 90.0) is None
    assert vcf._convert_het_to_hom(genos, info, "COV", 90.0) == 1
    assert vcf._convert_het_to_hom(genos, info, "COV", 90.1) is None
    info = {"COV": "90,9,1"}
    assert vcf._convert_het_to_hom(genos, info, "COV", 90.0) == 0
    assert vcf._convert_het_to_hom(genos, info, "COV", 90.1) is None


def test_bed_mask_file_to_dict():
    infile = os.path.join(data_dir, "bed_mask_file_to_dict.bed")
    expect = {"chrom1": [(0, 42)], "chrom2": [(20, 29), (100, 199)]}
    assert vcf._bed_mask_file_to_dict(infile) == expect


def test_vcf_to_variant_positions_to_mask_from_bed_file():
    vcf_file = os.path.join(
        data_dir, "vcf_to_variant_indexes_to_mask_from_bed_file.vcf"
    )
    bed_file = os.path.join(
        data_dir, "vcf_to_variant_indexes_to_mask_from_bed_file.bed"
    )
    expect = {"ref_42": {99}, "ref_43": {49, 50}, "ref_44": {18}}
    got = vcf.vcf_to_variant_positions_to_mask_from_bed_file(vcf_file, bed_file)
    assert got == expect


def test_load_vcf_file_for_distance_calc():
    infile = os.path.join(data_dir, "load_vcf_file_for_distance_calc.vcf")
    got_genos, got_counts = vcf.load_vcf_file_for_distance_calc(infile)
    expect_genos = np.array([0, 2, 2, 3, 4], dtype=np.uint16)
    np.testing.assert_array_equal(got_genos, expect_genos)
    expect_counts = variant_counts.VariantCounts(het=0, hom=3, null=1, het_to_hom=1)
    assert got_counts == expect_counts

    mask = {"ref_42": {10}}
    got_genos, got_counts = vcf.load_vcf_file_for_distance_calc(infile, mask=mask)
    expect_genos = np.array([2, 2, 3, 4], dtype=np.uint16)
    np.testing.assert_array_equal(got_genos, expect_genos)
    expect_counts = variant_counts.VariantCounts(het=0, hom=3, null=0, het_to_hom=1)
    assert got_counts == expect_counts

    got_genos, got_counts = vcf.load_vcf_file_for_distance_calc(
        infile, only_use_pass=False, het_to_hom_min_pc_depth=99.0
    )
    expect_genos = np.array([1, 0, 2, 3, 4], dtype=np.uint16)
    np.testing.assert_array_equal(got_genos, expect_genos)
    expect_counts = variant_counts.VariantCounts(het=1, hom=4, null=0, het_to_hom=0)
    assert got_counts == expect_counts

    got_genos, got_counts = vcf.load_vcf_file_for_distance_calc(
        infile,
        only_use_pass=False,
        numeric_filters={"GT_CONF": (True, 12)},
        het_to_hom_min_pc_depth=99.0,
    )
    expect_genos = np.array([1, 0, 0, 0, 4], dtype=np.uint16)
    np.testing.assert_array_equal(got_genos, expect_genos)
    expect_counts = variant_counts.VariantCounts(het=1, hom=2, null=2, het_to_hom=0)
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

    got = vcf.load_vcf_files_for_distance_calc(
        filenames, threads=2, only_use_pass=True, het_to_hom_key="ignore"
    )
    expect = [
        (
            np.array([0, 0, 2, 3, 4], dtype=np.uint16),
            variant_counts.VariantCounts(het=1, hom=3, null=1, het_to_hom=0),
        ),
        (
            np.array([0, 2, 2, 0, 2], dtype=np.uint16),
            variant_counts.VariantCounts(het=0, hom=3, null=2, het_to_hom=0),
        ),
    ]
    check_got_equal_expect(got, expect)

    got = vcf.load_vcf_files_for_distance_calc(
        filenames, threads=2, only_use_pass=False, het_to_hom_key="ignore"
    )
    expect = [
        (
            np.array([1, 0, 2, 3, 4], dtype=np.uint16),
            variant_counts.VariantCounts(het=1, hom=4, null=0, het_to_hom=0),
        ),
        (
            np.array([1, 2, 2, 3, 2], dtype=np.uint16),
            variant_counts.VariantCounts(het=0, hom=5, null=0, het_to_hom=0),
        ),
    ]
    check_got_equal_expect(got, expect)

    got = vcf.load_vcf_files_for_distance_calc(
        filenames,
        threads=2,
        only_use_pass=False,
        numeric_filters={"GT_CONF": (True, 12)},
        het_to_hom_key="ignore",
        mask_bed_file=os.path.join(
            data_dir, "load_vcf_files_for_distance_calc.mask.bed"
        ),
    )
    expect = [
        (
            np.array([0, 0, 0, 4], dtype=np.uint16),
            variant_counts.VariantCounts(het=1, hom=1, null=2, het_to_hom=0),
        ),
        (
            np.array([2, 2, 0, 2], dtype=np.uint16),
            variant_counts.VariantCounts(het=0, hom=3, null=1, het_to_hom=0),
        ),
    ]
    check_got_equal_expect(got, expect)


def test_sample_name_from_vcf():
    good_file = os.path.join(data_dir, "sample_name_from_vcf.good.vcf")
    assert vcf.sample_name_from_vcf(good_file) == "sample_42"
    bad_file = os.path.join(data_dir, "sample_name_from_vcf.bad.vcf")
    with pytest.raises(RuntimeError):
        vcf.sample_name_from_vcf(bad_file)


def test_sample_names_tsv_from_vcf_file_of_filenames():
    tmp_fofn = "tmp.test.sample_names_tsv_from_vcf_file_of_filenames.in"
    tmp_expect = "tmp.test.sample_names_tsv_from_vcf_file_of_filenames.expect"
    with open(tmp_fofn, "w") as f_fofn, open(tmp_expect, "w") as f_expect:
        print("sample\tvcf_file", file=f_expect)
        for i in range(1, 4, 1):
            filename = os.path.join(
                data_dir, f"sample_names_tsv_from_vcf_file_of_filenames.{i}.vcf"
            )
            print(filename, file=f_fofn)
            print(f"sample_{i}", os.path.abspath(filename), sep="\t", file=f_expect)

    tmp_out = "tmp.test.sample_names_tsv_from_vcf_file_of_filenames.out"
    if os.path.exists(tmp_out):
        os.unlink(tmp_out)
    vcf.sample_names_tsv_from_vcf_file_of_filenames(tmp_fofn, tmp_out, threads=2)
    assert filecmp.cmp(tmp_out, tmp_expect, shallow=False)
    os.unlink(tmp_fofn)
    os.unlink(tmp_expect)
    os.unlink(tmp_out)
