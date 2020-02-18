import os
import logging
import pytest

from triphecta import (
    genotypes,
    phenotypes,
    phenotype_compare,
    sample_neighbours_finding,
    strain_triple,
    vcf,
)

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "strain_triple")


@pytest.fixture(scope="function")
def genos(request):
    genos = genotypes.Genotypes(testing=True)
    genos.sample_names_list = ["s1", "s2", "s3", "s4", "s5", "s6", "s7"]
    genos._make_sample_name_to_index()
    genos.distances = {
        (0, 1): 0,
        (0, 2): 1,
        (0, 3): 1,
        (0, 4): 3,
        (0, 5): 5,
        (0, 6): 6,
        (1, 2): 6,
        (1, 3): 4,
        (1, 4): 7,
        (1, 5): 1,
        (1, 6): 2,
        (2, 3): 6,
        (2, 4): 4,
        (2, 5): 2,
        (2, 6): 3,
        (3, 4): 5,
        (3, 5): 8,
        (3, 6): 10,
        (4, 5): 4,
        (4, 6): 7,
        (5, 6): 1,
    }
    return genos


@pytest.fixture(scope="function")
def phenos(request):
    tmp_phenotypes_file = "tmp.phenotypes.tsv"
    with open(tmp_phenotypes_file, "w") as f:
        print("sample", "p1", "p2", "p3", "p4", "p5", "p6", sep="\t", file=f)
        print("s1", "R", "R", "R", "R", "S", "S", sep="\t", file=f)
        print("s2", "R", "R", "S", "S", "R", "R", sep="\t", file=f)
        print("s3", "S", "S", "R", "R", "S", "S", sep="\t", file=f)
        print("s4", "S", "S", "R", "R", "S", "S", sep="\t", file=f)
        print("s5", "S", "S", "R", "S", "R", "S", sep="\t", file=f)
        print("s6", "S", "S", "S", "S", "R", "R", sep="\t", file=f)
        print("s7", "S", "S", "S", "S", "R", "R", sep="\t", file=f)
    phenos = phenotypes.Phenotypes(tmp_phenotypes_file)
    os.unlink(tmp_phenotypes_file)
    return phenos


@pytest.fixture(scope="function")
def constraints(request):
    return {
        "p1": {"must_be_same": False, "method": "equal", "params": {}},
        "p2": {"must_be_same": False, "method": "equal", "params": {}},
        "p3": {"must_be_same": True, "method": "equal", "params": {}},
        "p4": {"must_be_same": True, "method": "equal", "params": {}},
        "p5": {"must_be_same": True, "method": "equal", "params": {}},
        "p6": {"must_be_same": True, "method": "equal", "params": {}},
    }


def test_find_strain_triples(genos, phenos, constraints, caplog):
    caplog.set_level(logging.INFO)
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    wanted_phenos = {"p1": "R", "p2": "R"}
    got = strain_triple.find_strain_triples(
        wanted_phenos, genos, phenos, pheno_compare, max_pheno_diffs=1, top_n_genos=10
    )
    rank_data_s1_1 = sample_neighbours_finding.RankData(
        sample="s3", rank_sum=0, geno_rank=0, pheno_rank=0, geno_dist=1, pheno_dist=0
    )
    rank_data_s1_2 = sample_neighbours_finding.RankData(
        sample="s4", rank_sum=0, geno_rank=0, pheno_rank=0, geno_dist=1, pheno_dist=0
    )
    rank_data_s2_1 = sample_neighbours_finding.RankData(
        sample="s6", rank_sum=0, geno_rank=0, pheno_rank=0, geno_dist=1, pheno_dist=0
    )
    rank_data_s2_2 = sample_neighbours_finding.RankData(
        sample="s7", rank_sum=1, geno_rank=1, pheno_rank=0, geno_dist=2, pheno_dist=0
    )
    expect = [
        strain_triple.StrainTriple("s1", rank_data_s1_1, rank_data_s1_2),
        strain_triple.StrainTriple("s2", rank_data_s2_1, rank_data_s2_2),
    ]
    assert got == expect


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
    expect_variant_calls = {"case": [0, 1], "control1": [0, 2], "control2": [0, 2]}
    assert triple.variant_calls == expect_variant_calls


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
        "case": [0, 0, ".", 0, 0],
        "control1": [1, ".", 1, 0, 1],
        "control2": [".", 1, 1, 1, 1],
    }
    assert triple.variant_indexes_of_interest == set()
    triple.update_variants_of_interest()
    assert triple.variant_indexes_of_interest == {4}
