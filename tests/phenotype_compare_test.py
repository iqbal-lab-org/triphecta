import pytest

from triphecta import phenotype_compare


def test_init_in_particular_sanity_check_constraints():
    constraints = {"d1": {"must_be_same": True, "method": "equal", "params": {}}}
    phenotype_compare.PhenotypeCompare(constraints)

    # Fails requirement that must be at least one where "muste_be_same" is True
    constraints = {"d1": {"must_be_same": False, "method": "equal", "params": {}}}
    with pytest.raises(RuntimeError):
        phenotype_compare.PhenotypeCompare(constraints)

    # Uses an unknown method
    constraints = {"d1": {"must_be_same": True, "method": "WRONG", "params": {}}}
    with pytest.raises(RuntimeError):
        phenotype_compare.PhenotypeCompare(constraints)

    # equal method should not have any params
    constraints = {
        "d1": {"must_be_same": True, "method": "equal", "params": {"low": 1}}
    }
    with pytest.raises(RuntimeError):
        phenotype_compare.PhenotypeCompare(constraints)

    # Range method should have low and high params
    constraints = {"d1": {"must_be_same": True, "method": "range", "params": {}}}
    with pytest.raises(RuntimeError):
        phenotype_compare.PhenotypeCompare(constraints)
    constraints = {
        "d1": {"must_be_same": True, "method": "range", "params": {"low": 1}}
    }
    with pytest.raises(RuntimeError):
        phenotype_compare.PhenotypeCompare(constraints)
    constraints = {
        "d1": {"must_be_same": True, "method": "range", "params": {"low": 1, "high": 2}}
    }
    phenotype_compare.PhenotypeCompare(constraints)


def test_satisfy_required_differences_one_the_same():
    constraints = {
        "d0": {"must_be_same": True, "method": "equal", "params": {}},
        "d1": {"must_be_same": False, "method": "equal", "params": {}},
    }
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    p1 = {"d1": True}
    p2 = {"d1": True}
    p3 = {"d1": False}
    p4 = {"d1": None}
    p5 = {"d1": None}

    assert pheno_compare.satisfy_required_differences(p1, p3)
    assert pheno_compare.satisfy_required_differences(p3, p1)
    assert not pheno_compare.satisfy_required_differences(p1, p1)
    assert not pheno_compare.satisfy_required_differences(p1, p2)
    assert not pheno_compare.satisfy_required_differences(p2, p1)
    assert not pheno_compare.satisfy_required_differences(p1, p4)
    assert not pheno_compare.satisfy_required_differences(p4, p1)
    assert not pheno_compare.satisfy_required_differences(p4, p5)
    assert not pheno_compare.satisfy_required_differences(p5, p4)


def test_satisfy_required_differences_two_the_same():
    constraints = {
        "d0": {"must_be_same": True, "method": "equal", "params": {}},
        "d1": {"must_be_same": True, "method": "equal", "params": {}},
        "d2": {
            "must_be_same": False,
            "method": "range",
            "params": {"low": 1, "high": 2},
        },
    }
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    p1 = {"d0": True, "d1": True, "d2": 1.5}
    p2 = {"d0": True, "d1": True, "d2": 1.1}
    p3 = {"d0": True, "d1": False, "d2": 1.5}
    p4 = {"d0": True, "d1": True, "d2": 1.5}
    p5 = {"d0": True, "d1": True, "d2": 2.5}
    assert not pheno_compare.satisfy_required_differences(p1, p2)
    assert not pheno_compare.satisfy_required_differences(p2, p1)
    assert not pheno_compare.satisfy_required_differences(p1, p3)
    assert not pheno_compare.satisfy_required_differences(p3, p1)
    assert not pheno_compare.satisfy_required_differences(p1, p4)
    assert not pheno_compare.satisfy_required_differences(p4, p1)
    assert pheno_compare.satisfy_required_differences(p1, p5)
    assert pheno_compare.satisfy_required_differences(p5, p1)


def test_phenos_equal_account_for_none_using_method_equal():
    f = phenotype_compare.PhenotypeCompare._phenos_equal_account_for_none
    c = phenotype_compare.PhenotypeCompare._compare_method_equal
    assert f(True, True, c, True)
    assert f(True, True, c, False)
    assert f(False, False, c, True)
    assert f(False, False, c, False)
    assert not f(True, False, c, True)
    assert not f(True, False, c, False)
    assert not f(False, True, c, True)
    assert not f(False, True, c, False)
    assert f("x", None, c, False)
    assert not f("x", None, c, True)
    assert f("R", "R", c, True)
    assert f("R", "R", c, False)
    assert not f("R", "S", c, True)
    assert not f("R", "S", c, False)
    assert f("R", None, c, False)
    assert not f("R", None, c, True)
    assert not f(None, None, c, True)
    assert f(None, None, c, False)


def test_phenos_equal_account_for_none_using_method_range():
    f = phenotype_compare.PhenotypeCompare._phenos_equal_account_for_none
    c = phenotype_compare.PhenotypeCompare._compare_method_range
    assert f(1, 1.1, c, True, low=1, high=2)
    assert f(1, 1.1, c, False, low=1, high=2)
    assert not f(1, 1.5, c, True, low=1.25, high=2)
    assert not f(1, 1.5, c, False, low=1.25, high=2)
    assert not f(1.5, 2.5, c, True, low=1.25, high=2)
    assert not f(1.5, 2.5, c, False, low=1.25, high=2)
    assert f(0, 3, c, True, low=1, high=2)
    assert f(0, 3, c, False, low=1, high=2)
    assert not f(None, 1.1, c, True, low=1, high=2)
    assert f(None, 1.1, c, False, low=1, high=2)
    assert not f(None, 0, c, True, low=1, high=2)
    assert f(None, 0, c, False, low=1, high=2)


def test_differences():
    constraints = {
        "d0": {"must_be_same": True, "method": "equal", "params": {}},
        "d1": {"must_be_same": True, "method": "equal", "params": {}},
        "d2": {"must_be_same": True, "method": "equal", "params": {}},
        "d3": {
            "must_be_same": False,
            "method": "range",
            "params": {"low": 1, "high": 2},
        },
    }
    pheno_compare = phenotype_compare.PhenotypeCompare(constraints)
    p1 = {"d0": True, "d1": True, "d2": True, "d3": 2.5}
    p2 = {"d0": True, "d1": True, "d2": True, "d3": 2.5}
    p3 = {"d0": True, "d1": True, "d2": False, "d3": 2.5}
    p4 = {"d0": True, "d1": False, "d2": False, "d3": 2.5}
    p5 = {"d0": False, "d1": False, "d2": False, "d3": 2.5}
    assert pheno_compare.differences(p1, p1) == 0
    assert pheno_compare.differences(p1, p2) == 0
    assert pheno_compare.differences(p2, p1) == 0
    assert pheno_compare.differences(p1, p3) == 1
    assert pheno_compare.differences(p3, p1) == 1
    assert pheno_compare.differences(p1, p4) == 2
    assert pheno_compare.differences(p4, p1) == 2
    assert pheno_compare.differences(p1, p5) == 3
    assert pheno_compare.differences(p5, p1) == 3
