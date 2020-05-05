import collections
import filecmp
import os
import logging
import pytest
import subprocess

from triphecta import (
    genotypes,
    phenotypes,
    phenotype_compare,
    sample_neighbours_finding,
    strain_triple,
    strain_triples,
    vcf,
)

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "strain_triples")


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
    genos.vcf_files = {
        "s1": os.path.join(data_dir, "genos.vcf.1"),
        "s2": os.path.join(data_dir, "genos.vcf.1"),
        "s3": os.path.join(data_dir, "genos.vcf.1"),
        "s4": os.path.join(data_dir, "genos.vcf.2"),
        "s5": os.path.join(data_dir, "genos.vcf.2"),
        "s6": os.path.join(data_dir, "genos.vcf.2"),
        "s7": os.path.join(data_dir, "genos.vcf.3"),
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
    triples = strain_triples.StrainTriples(genos, phenos, pheno_compare, top_n_genos=10)
    case_sample_names = ["s1", "s2"]
    triples.find_strain_triples(case_sample_names)

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

    assert triples.triples == expect


def test_write_triples_names_file():
    tmp_out = "tmp.strain_triples.write_triples_names_file.tsv"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    Triple = collections.namedtuple("Triple", ["case", "control1", "control2"])
    Control = collections.namedtuple("Control", ["sample", "geno_dist", "pheno_dist"])
    triples = [
        Triple("case1", Control("control1.1", 1, 2), Control("control1.2", 4, 5)),
        Triple("case2", Control("control2.1", 2, 3), Control("control2.2", 5, 6)),
        Triple("case3", Control("control3.2", 3, 4), Control("control3.2", 6, 7)),
    ]
    strain_triples.StrainTriples._write_triples_names_file(triples, tmp_out)
    expect = os.path.join(data_dir, "write_triples_names_file.tsv")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)


def test_write_variants_summary_file():
    triples = [
        strain_triple.StrainTriple("case1", "control1", "control2"),
        strain_triple.StrainTriple("case2", "control2", "control2"),
        strain_triple.StrainTriple("case3", "control3", "control3"),
        strain_triple.StrainTriple("case4", "control4", "control4"),
    ]
    triples[0].variants = [
        vcf.Variant(CHROM="ref_1", POS=9, REF="A", ALTS=["C"]),
        vcf.Variant(CHROM="ref_1", POS=10, REF="A", ALTS=["C,G"]),
        vcf.Variant(CHROM="ref_1", POS=11, REF="C", ALTS=["T"]),
        vcf.Variant(CHROM="ref_2", POS=42, REF="T", ALTS=["A"]),
    ]
    triples[0].variant_indexes_of_interest = {0, 1, 2, 3}
    triples[1].variant_indexes_of_interest = set()
    triples[2].variant_indexes_of_interest = {0, 1, 2}
    triples[3].variant_indexes_of_interest = {2}
    tmp_out = "tmp.strain_triples.write_variants_summary_file.tsv"
    subprocess.check_output(f"rm -f {tmp_out}", shell=True)
    strain_triples.StrainTriples._write_variants_summary_file(
        triples, tmp_out, vcf_records_to_mask=None
    )
    expect = os.path.join(data_dir, "write_variants_summary_file.no_mask.tsv")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)

    mask = {"ref_1": {9, 11}, "ref_2": {42}}
    strain_triples.StrainTriples._write_variants_summary_file(
        triples, tmp_out, vcf_records_to_mask=mask
    )
    expect = os.path.join(data_dir, "write_variants_summary_file.with_mask.tsv")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)


def test_run_analysis(genos, phenos, constraints, caplog):
    caplog.set_level(logging.INFO)
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    triples = strain_triples.StrainTriples(genos, phenos, pheno_compare, top_n_genos=10)
    case_sample_names = ["s1", "s2"]
    tmp_dir = "tmp.strain_triples.run_analysis"
    subprocess.check_output(f"rm -rf {tmp_dir}", shell=True)
    os.mkdir(tmp_dir)
    outprefix = os.path.join(tmp_dir, "out")
    mask_file = os.path.join(data_dir, "run_analysis.mask.bed")
    got = triples.run_analysis(case_sample_names, outprefix, mask_file=mask_file)
    # The contents of these files are checked elsewhere. This test is just to
    # check then run of whole pipeline doesn't crash
    for filename in got.values():
        assert os.path.exists(filename)
    subprocess.check_output(f"rm -r {tmp_dir}", shell=True)
