import os

import pytest

from triphecta import distances

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "distances")


def test_distances_between_vcf_files():
    filenames = {
        "s1": os.path.join(data_dir, "distances_between_vcf_files.1.vcf"),
        "s2": os.path.join(data_dir, "distances_between_vcf_files.2.vcf"),
        "s3": os.path.join(data_dir, "distances_between_vcf_files.3.vcf"),
    }

    got_sample_names, got_dists, got_variant_counts = distances.distances_between_vcf_files(
        filenames, threads=2
    )
    expect_sample_names = ["s1", "s2", "s3"]
    expect_dists = {(0, 1): 1, (0, 2): 2, (1, 2): 0}
    expect_variant_counts = [
        {"hom": 3, "het": 1, "null": 1},
        {"hom": 3, "het": 0, "null": 2},
        {"hom": 3, "het": 1, "null": 1},
    ]

    assert got_sample_names == expect_sample_names
    assert got_dists == expect_dists
    assert got_variant_counts == expect_variant_counts

    pickle_file = "tmp.test_distances_between_vcf_files.pickle"
    if os.path.exists(pickle_file):
        os.unlink(pickle_file)

    file_of_filenames = "tmp.distances_between_vcf_files.fofn"
    with open(file_of_filenames, "w") as f:
        print("sample\tvcf_file", file=f)
        for sample, filename in filenames.items():
            print(sample, filename, sep="\t", file=f)
    distances.pickle_distances_between_vcf_files(
        file_of_filenames, pickle_file, threads=2
    )
    os.unlink(file_of_filenames)
    got_sample_names, got_dists, got_variant_counts = distances.load_from_pickle(pickle_file)
    assert got_dists == expect_dists
    assert got_variant_counts == expect_variant_counts
    if os.path.exists(pickle_file):
        os.unlink(pickle_file)


def test_load_one_sample_distances_file():
    dist_file = os.path.join(data_dir, "load_one_sample_distances_file.tsv")
    expect = [("s1", 0.0), ("s2", 42.0), ("s3", 100.0)]
    got = distances._load_one_sample_distances_file(dist_file)
    assert got == expect


def test_update_distances_for_one_sample():
    sample_names = ["s1", "s2", "s3"]
    sample_names_to_index = {sample_names[i]: i for i in range(len(sample_names))}
    data1 = [("s1", 0.0), ("s2", 42.0), ("s3", 100.0)]
    data2 = [("s1", 42.0), ("s3", 50.0)]
    data3 = [("s1", 200.0)]
    all_dists = {}
    distances._update_distances_for_one_sample(0, data1, all_dists, sample_names_to_index)
    expect_distances = {(0, 1): 42, (0, 2): 100}
    assert expect_distances == all_dists

    distances._update_distances_for_one_sample(1, data2, all_dists, sample_names_to_index)
    expect_distances[(1, 2)] = 50
    assert expect_distances == all_dists

    # This has distance s3 to s1 of 200, which doesn't agree with the
    # existing value of 100
    with pytest.raises(RuntimeError):
        distances._update_distances_for_one_sample(2, data3, all_dists, sample_names_to_index)


def test_load_sample_distances_file_of_filenames():
    infile = os.path.join(data_dir, "load_sample_distances_file_of_filenames.tsv")
    expect_names = ["s1", "s2"]
    expect_files = ["f1", "f2"]
    got_names, got_files = distances._load_sample_distances_file_of_filenames(infile)
    assert got_names == ["s1", "s2"]
    assert got_files == ["f1", "f2"]


def test_load_all_one_sample_distances_files():
    tmp_tsv = "tmp.test_load_all_one_sample_distances_files.tsv"
    with open(tmp_tsv, "w") as f:
        print("sample", "distance_file", sep="\t", file=f)
        for i in range(4):
            print(
                f"s{i+1}",
                os.path.join(
                    data_dir, f"load_all_one_sample_distances_files.{i+1}.tsv"
                ),
                sep="\t",
                file=f,
            )
    got_names, got_dists = distances.load_all_one_sample_distances_files(tmp_tsv, threads=2)
    expect_names = ["s1", "s2", "s3", "s4"]
    expect_dists = {
        (0, 1): 3.0,
        (0, 2): 4.0,
        (0, 3): 6.0,
        (1, 2): 10.0,
        (1, 3): 8.0,
        (2, 3): 5.0,
    }
    assert got_names == expect_names
    assert got_dists == expect_dists

    pickle_file = "tmp.test_load_all_one_sample_distances_files.pickle"
    if os.path.exists(pickle_file):
        os.unlink(pickle_file)
    distances.pickle_load_all_one_sample_distances_files(
        tmp_tsv, pickle_file, threads=2
    )
    got_names, got_dists, got_variant_counts = distances.load_from_pickle(pickle_file)
    assert got_dists == expect_dists
    assert got_variant_counts == {}
    if os.path.exists(pickle_file):
        os.unlink(pickle_file)
    os.unlink(tmp_tsv)
