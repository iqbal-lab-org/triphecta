import filecmp
import os
import logging
import pytest

from triphecta import strain_triple, vcf

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "strain_triple")


def test_load_variants_from_vcf_files():
    triple = strain_triple.StrainTriple("case", "control1", "control2")
    case_vcf = os.path.join(data_dir, "load_variants_from_vcf_files.case.vcf")
    control1_vcf = os.path.join(data_dir, "load_variants_from_vcf_files.control1.vcf")
    control2_vcf = os.path.join(data_dir, "load_variants_from_vcf_files.control2.vcf")
    triple.load_variants_from_vcf_files(case_vcf, control1_vcf, control2_vcf)
    expect_variants = [
        vcf.Variant(CHROM="ref_42", POS=10, REF="C", ALTS=["G"]),
        vcf.Variant(CHROM="ref_43", POS=41, REF="T", ALTS=["A", "CT"]),
    ]
    assert triple.variants == expect_variants
    expect_variant_calls = {"case": [{0}, {1}], "control1": [{0}, {2}], "control2": [{0}, {2}]}
    assert triple.variant_calls == expect_variant_calls

    triple = strain_triple.StrainTriple("case", "control1", "control2")
    triple.set_variants(expect_variants)
    triple.load_variants_from_vcf_files(case_vcf, control1_vcf, control2_vcf)
    assert triple.variants == expect_variants
    assert triple.variant_calls == expect_variant_calls
    triple.clear_variant_calls()
    assert triple.variant_calls == {"case": None, "control1": None, "control2": None}



def test_genotypes_are_of_interest():
    f = strain_triple.StrainTriple.genotypes_are_of_interest
    assert not f(None, None, None)
    assert not f({0}, {1}, None)
    assert not f({0}, None, {1})
    assert f({0}, {1}, {1})
    assert f({1}, {0}, {0})
    assert not f({0, 1}, {2, 3}, {4, 5})
    assert not f({0, 1}, {2}, {3})
    assert f({0, 1}, {2}, {2})
    assert f({0, 1}, {2}, {2, 3})
    assert f({1}, {2}, {2, 3})
    assert not f({1}, {2}, {1, 2})
    assert f({1}, {2, 3}, {3, 4})
    assert f({0, 1}, {2, 3}, {3, 4})
    assert not f({0, 2}, {2, 3}, {3, 4})
    assert not f({0, 3}, {2, 3}, {3, 4})
    assert not f({0, 4}, {2, 3}, {3, 4})
    assert not f({0, 1}, {0, 1}, {0, 1})
    assert f({0}, {3, 4}, {3, 4})
    assert f({0, 1}, {3, 4}, {3, 4})


def test_update_variants_of_interest():
    triple = strain_triple.StrainTriple("case", "control1", "control2")
    triple.variants = [
        vcf.Variant(CHROM="ref_42", POS=9, REF="A", ALTS=["C"]),
        vcf.Variant(CHROM="ref_42", POS=10, REF="A", ALTS=["C,G"]),
        vcf.Variant(CHROM="ref_42", POS=11, REF="C", ALTS=["A,G"]),
        vcf.Variant(CHROM="ref_42", POS=12, REF="G", ALTS=["T"]),
        vcf.Variant(CHROM="ref_42", POS=13, REF="T", ALTS=["A"]),
    ]
    triple.variant_calls = {
        "case": [{0}, {0}, None, {0}, {0}],
        "control1": [{1}, None, {1}, {0}, {1}],
        "control2": [None, {1}, {1}, {1}, {1}],
    }
    assert triple.variant_indexes_of_interest == set()
    triple.update_variants_of_interest()
    assert triple.variant_indexes_of_interest == {4}


def test_genotype_to_string():
    f = strain_triple.StrainTriple.genotype_to_string
    assert f(None) == "./."
    assert f({0}) == "0/0"
    assert f({0, 1}) == "0/1"


def test_write_variants_of_interest_file():
    triple = strain_triple.StrainTriple("case", "control1", "control2")
    triple.variants = [
        vcf.Variant(CHROM="ref_42", POS=9, REF="A", ALTS=["C"]),
        vcf.Variant(CHROM="ref_42", POS=10, REF="A", ALTS=["C,G"]),
        vcf.Variant(CHROM="ref_42", POS=11, REF="C", ALTS=["A,G"]),
        vcf.Variant(CHROM="ref_42", POS=12, REF="G", ALTS=["T"]),
        vcf.Variant(CHROM="ref_42", POS=13, REF="T", ALTS=["A"]),
        vcf.Variant(CHROM="ref_43", POS=42, REF="A", ALTS=["C,G"]),
    ]
    triple.variant_calls = {
        "case": [{0}, {0}, None, {0}, {0}, {0, 1}],
        "control1": [{1}, None, {1}, {0}, {1}, {2}],
        "control2": [None, {1}, {1}, {1}, {1}, {2}],
    }
    triple.variant_indexes_of_interest = {4,5}

    outfile = "tmp.out.write_variants_of_interest_file.tsv"
    if os.path.exists(outfile):
        os.unlink(outfile)
    triple.write_variants_of_interest_file(outfile)
    expect = os.path.join(data_dir, "write_variants_of_interest_file.tsv")
    assert filecmp.cmp(expect, outfile, shallow=False)

    mask = {"ref_42": {13}}
    triple.write_variants_of_interest_file(outfile, vcf_records_to_mask=mask)
    expect = os.path.join(data_dir, "write_variants_of_interest_file.masked.tsv")
    assert filecmp.cmp(expect, outfile, shallow=False)
    os.unlink(outfile)

