import filecmp
import os
import pytest

from triphecta import tree, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "tree")


def test_newick_from_dist_matrix():
    infile = os.path.join(data_dir, "newick_from_dist_matrix.txt")
    tmp_out = "tmp.newick_from_dist_matrix.out"
    utils.rm_rf(tmp_out)
    tree.newick_from_dist_matrix(infile, tmp_out, "upgma")
    assert os.path.exists(tmp_out)
    os.unlink(tmp_out)
    tree.newick_from_dist_matrix(infile, tmp_out, "nj")
    assert os.path.exists(tmp_out)
    os.unlink(tmp_out)
