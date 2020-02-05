import os

import pytest

from triphecta import genotypes

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genotypes")


def test_distance():
    genos = genotypes.Genotypes(testing=True)
    genos.vcf_files = {"s1", "s2", "s3"}
    genos.distances = {("s1", "s2"): 0, ("s1", "s3"): 3, ("s2", "s3"): 5}
    assert genos.distance("s1", "s2") == 0
    assert genos.distance("s2", "s1") == 0
    assert genos.distance("s1", "s3") == 3
    assert genos.distance("s3", "s1") == 3
    assert genos.distance("s2", "s3") == 5
    assert genos.distance("s3", "s2") == 5


def test_distance_dict():
    genos = genotypes.Genotypes(testing=True)
    genos.vcf_files = {"s1", "s2", "s3", "s4", "s5", "s6", "s7"}
    genos.distances = {
        ("s1", "s2"): 0,
        ("s1", "s3"): 3,
        ("s1", "s4"): 3,
        ("s1", "s5"): 4,
        ("s1", "s6"): 4,
        ("s1", "s7"): 5,
    }
    assert genos.distance_dict("s1") == {
        "s2": 0,
        "s3": 3,
        "s4": 3,
        "s5": 4,
        "s6": 4,
        "s7": 5,
    }
    assert genos.distance_dict("s1", top_n=1) == {"s2": 0}
    assert genos.distance_dict("s1", top_n=2) == {"s2": 0, "s3": 3, "s4": 3}
    assert genos.distance_dict("s1", top_n=3) == {"s2": 0, "s3": 3, "s4": 3}
    assert genos.distance_dict("s1", top_n=4) == {
        "s2": 0,
        "s3": 3,
        "s4": 3,
        "s5": 4,
        "s6": 4,
    }
    assert genos.distance_dict("s1", top_n=5) == {
        "s2": 0,
        "s3": 3,
        "s4": 3,
        "s5": 4,
        "s6": 4,
    }
    assert genos.distance_dict("s1", top_n=6) == {
        "s2": 0,
        "s3": 3,
        "s4": 3,
        "s5": 4,
        "s6": 4,
        "s7": 5,
    }
