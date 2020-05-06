import os
import pytest

from triphecta import genotypes, variant_counts

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genotypes")


def test_distance():
    genos = genotypes.Genotypes(testing=True)
    genos.sample_names_list = ["s1", "s2", "s3"]
    genos._make_sample_name_to_index()
    genos.distances = {(0, 1): 0, (0, 2): 3, (1, 2): 5}
    assert genos.distance("s1", "s2") == 0
    assert genos.distance("s2", "s1") == 0
    assert genos.distance("s1", "s3") == 3
    assert genos.distance("s3", "s1") == 3
    assert genos.distance("s2", "s3") == 5
    assert genos.distance("s3", "s2") == 5


def test_distance_dict():
    genos = genotypes.Genotypes(testing=True)
    genos.sample_names_list = ["s1", "s2", "s3", "s4", "s5", "s6", "s7"]
    genos._make_sample_name_to_index()
    genos.distances = {(0, 1): 0, (0, 2): 3, (0, 3): 3, (0, 4): 4, (0, 5): 4, (0, 6): 5}
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

    genos.excluded_samples = {"s2": "reason"}
    assert genos.distance_dict("s1", top_n=1) == {"s3": 3, "s4": 3}


def test_update_excluded_samples_using_variant_counts():
    genos = genotypes.Genotypes(testing=True)
    assert genos.excluded_samples == {}
    genos.sample_names_list = ["s1", "s2", "s3", "s4"]
    genos._make_sample_name_to_index()
    genos.vcf_variant_counts = [
        variant_counts.VariantCounts(hom=95, het=5, null=0, het_to_hom=0),
        variant_counts.VariantCounts(hom=90, het=10, null=0, het_to_hom=0),
        variant_counts.VariantCounts(hom=89, het=11, null=0, het_to_hom=0),
        variant_counts.VariantCounts(hom=89, het=10, null=0, het_to_hom=1),
        variant_counts.VariantCounts(hom=88, het=10, null=1, het_to_hom=1),
    ]
    genos.update_excluded_samples_using_variant_counts()
    expect = {x: {"Too few hom calls"} for x in [2, 4]}
    assert genos.excluded_samples == expect
    genos.excluded_samples = {}
    genos.update_excluded_samples_using_variant_counts(minimum_percent_hom_calls=95)
    expect = {x: {"Too few hom calls"} for x in [1, 2, 3, 4]}
    assert genos.excluded_samples == expect
