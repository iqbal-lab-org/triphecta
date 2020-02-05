import os

import pytest

from triphecta import (
    genotypes,
    phenotypes,
    phenotype_compare,
    sample_neighbours_finding,
)

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "sample_neighbours_finding")


@pytest.fixture(scope="function")
def genos(request):
    genos = genotypes.Genotypes(testing=True)
    genos.vcf_files = {"s1", "s2", "s3", "s4", "s5"}
    genos.distances = {
        ("s1", "s2"): 0,
        ("s1", "s3"): 3,
        ("s1", "s4"): 2,
        ("s1", "s5"): 1,
        ("s2", "s3"): 1,
        ("s2", "s4"): 4,
        ("s2", "s5"): 3,
        ("s3", "s4"): 5,
        ("s3", "s5"): 2,
        ("s4", "s4"): 5,
        ("s4", "s5"): 1,
    }
    return genos


@pytest.fixture(scope="function")
def phenos(request):
    tmp_phenotypes_file = "tmp.phenotypes.tsv"
    with open(tmp_phenotypes_file, "w") as f:
        print("sample", "p1", "p2", "p3", "p4", sep="\t", file=f)
        print("s1", "T", "T", "T", "T", sep="\t", file=f)
        print("s2", "T", "T", "T", "F", sep="\t", file=f)
        print("s3", "F", "T", "T", "F", sep="\t", file=f)
        print("s4", "F", "T", "F", "F", sep="\t", file=f)
        print("s5", "F", "F", "F", "F", sep="\t", file=f)
    phenos = phenotypes.Phenotypes(tmp_phenotypes_file)
    os.unlink(tmp_phenotypes_file)
    return phenos


@pytest.fixture(scope="function")
def constraints(request):
    return {
        "p1": {"must_be_same": False, "method": "equal", "params": {}},
        "p2": {"must_be_same": True, "method": "equal", "params": {}},
        "p3": {"must_be_same": True, "method": "equal", "params": {}},
        "p4": {"must_be_same": True, "method": "equal", "params": {}},
    }


def test_geno_and_pheno_distances_for_one_sample(genos, phenos, constraints):
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)

    got_geno_distances, got_pheno_distances = sample_neighbours_finding._geno_and_pheno_distances_for_one_sample(
        genos, phenos, pheno_compare, "s1"
    )
    expect_geno_distances = {"s3": 3, "s4": 2, "s5": 1}
    expect_pheno_distances = {"s3": 1, "s4": 2, "s5": 3}
    assert got_geno_distances == expect_geno_distances
    assert got_pheno_distances == expect_pheno_distances

    got_geno_distances, got_pheno_distances = sample_neighbours_finding._geno_and_pheno_distances_for_one_sample(
        genos, phenos, pheno_compare, "s1", top_n_genos=3
    )
    expect_geno_distances = {"s4": 2, "s5": 1}
    expect_pheno_distances = {"s4": 2, "s5": 3}
    assert got_geno_distances == expect_geno_distances
    assert got_pheno_distances == expect_pheno_distances


def test_geno_and_pheno_distances_to_rank_table():
    geno_distances = {"s1": 1, "s2": 1, "s3": 5, "s4": 4}
    pheno_distances = {"s1": 0, "s2": 1, "s3": 2, "s4": 3}
    got = sample_neighbours_finding._geno_and_pheno_distances_to_rank_table(
        geno_distances, pheno_distances
    )
    expect = [
        sample_neighbours_finding.RankData(
            sample="s1",
            rank_sum=0,
            geno_rank=0,
            pheno_rank=0,
            geno_dist=1,
            pheno_dist=0,
        ),
        sample_neighbours_finding.RankData(
            sample="s2",
            rank_sum=1,
            geno_rank=0,
            pheno_rank=1,
            geno_dist=1,
            pheno_dist=1,
        ),
        sample_neighbours_finding.RankData(
            sample="s4",
            rank_sum=4,
            geno_rank=1,
            pheno_rank=3,
            geno_dist=4,
            pheno_dist=3,
        ),
        sample_neighbours_finding.RankData(
            sample="s3",
            rank_sum=4,
            geno_rank=2,
            pheno_rank=2,
            geno_dist=5,
            pheno_dist=2,
        ),
    ]
    assert expect == got


def test_ranked_neighbours_for_one_sample(genos, phenos, constraints):
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    got = sample_neighbours_finding.ranked_neighbours_for_one_sample(
        genos, phenos, pheno_compare, "s1"
    )
    expect = [
        sample_neighbours_finding.RankData(
            sample="s5",
            rank_sum=2,
            geno_rank=0,
            pheno_rank=2,
            geno_dist=1,
            pheno_dist=3,
        ),
        sample_neighbours_finding.RankData(
            sample="s4",
            rank_sum=2,
            geno_rank=1,
            pheno_rank=1,
            geno_dist=2,
            pheno_dist=2,
        ),
        sample_neighbours_finding.RankData(
            sample="s3",
            rank_sum=2,
            geno_rank=2,
            pheno_rank=0,
            geno_dist=3,
            pheno_dist=1,
        ),
    ]
    assert got == expect
    got = sample_neighbours_finding.ranked_neighbours_for_one_sample(
        genos, phenos, pheno_compare, "s1", top_n_genos=2
    )
    expect = [
        sample_neighbours_finding.RankData(
            sample="s5",
            rank_sum=0,
            geno_rank=0,
            pheno_rank=0,
            geno_dist=1,
            pheno_dist=3,
        )
    ]
    assert got == expect
