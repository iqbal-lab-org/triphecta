import filecmp
import os
import pytest

from triphecta import distances, utils, variant_counts, vcf

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "distances")


def test_distances_between_vcf_files():
    vcf_names_tsv = os.path.join(data_dir, "distances_between_vcf_files.vcfs.tsv")
    mask_bed_file = os.path.join(data_dir, "distances_between_vcf_files.mask.bed")
    outprefix = "tmp.distances_between_vcf_files"
    got_matrix_file = f"{outprefix}.distance_matrix.txt.gz"
    got_variant_counts_file = f"{outprefix}.variant_counts.tsv.gz"
    utils.rm_rf(got_matrix_file, got_variant_counts_file)

    (
        got_sample_names,
        got_dists,
        got_variant_counts,
    ) = distances.distances_between_vcf_files(
        vcf_names_tsv,
        outprefix,
        threads=2,
        het_to_hom_key="ignore",
        mask_bed_file=mask_bed_file,
    )
    expect_sample_names = ["s1", "s2", "s3"]
    expect_dists = {(0, 1): 1, (0, 2): 2, (1, 2): 0}
    expect_variant_counts = [
        variant_counts.VariantCounts(hom=3, het=1, null=1, het_to_hom=0),
        variant_counts.VariantCounts(hom=3, het=0, null=2, het_to_hom=0),
        variant_counts.VariantCounts(hom=3, het=1, null=1, het_to_hom=0),
    ]

    assert got_sample_names == expect_sample_names
    assert got_dists == expect_dists
    assert got_variant_counts == expect_variant_counts

    loaded_sample_names, loaded_dists = distances.load_distance_matrix_file(
        got_matrix_file
    )
    assert got_sample_names == loaded_sample_names
    assert got_dists == loaded_dists
    loaded_variant_counts = variant_counts.load_variant_count_list_from_tsv(
        got_variant_counts_file
    )
    assert got_variant_counts == loaded_variant_counts

    os.unlink(got_matrix_file)
    os.unlink(got_variant_counts_file)


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
    distances._update_distances_for_one_sample(
        0, data1, all_dists, sample_names_to_index
    )
    expect_distances = {(0, 1): 42, (0, 2): 100}
    assert expect_distances == all_dists

    distances._update_distances_for_one_sample(
        1, data2, all_dists, sample_names_to_index
    )
    expect_distances[(1, 2)] = 50
    assert expect_distances == all_dists

    # This has distance s3 to s1 of 200, which doesn't agree with the
    # existing value of 100
    with pytest.raises(RuntimeError):
        distances._update_distances_for_one_sample(
            2, data3, all_dists, sample_names_to_index
        )


def test_load_sample_distances_file_of_filenames():
    infile = os.path.join(data_dir, "load_sample_distances_file_of_filenames.tsv")
    expect_names = ["s1", "s2"]
    expect_files = ["f1", "f2"]
    got_names, got_files = distances._load_sample_distances_file_of_filenames(infile)
    assert got_names == ["s1", "s2"]
    assert got_files == ["f1", "f2"]


def test_distances_from_all_one_sample_distances_files():
    tmp_tsv = "tmp.test_distances_from_all_one_sample_distances_files.in.tsv"
    tmp_out = "tmp.test_distances_from_all_one_sample_distances_files.out.tsv"
    utils.rm_rf(tmp_out)
    with open(tmp_tsv, "w") as f:
        print("sample", "distance_file", sep="\t", file=f)
        for i in range(4):
            print(
                f"s{i+1}",
                os.path.join(
                    data_dir, f"distances_from_all_one_sample_distances_files.{i+1}.tsv"
                ),
                sep="\t",
                file=f,
            )
    got_names, got_dists = distances.distances_from_all_one_sample_distances_files(
        tmp_tsv, tmp_out, threads=2
    )
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

    loaded_names, loaded_dists = distances.load_distance_matrix_file(tmp_out)
    assert loaded_names == expect_names
    assert loaded_dists == expect_dists
    os.unlink(tmp_tsv)
    os.unlink(tmp_out)


def test_write_distance_matrix_file():
    tmp_out = "tmp.distances.write_distance_matrix_file.txt"
    utils.rm_rf(tmp_out)
    sample_names = ["sample1", "sample2", "sample3"]
    dists = {(0, 1): 3, (0, 2): 4, (1, 2): 42}
    distances.write_distance_matrix_file(sample_names, dists, tmp_out)
    expect = os.path.join(data_dir, "write_distance_matrix_file.txt")
    assert filecmp.cmp(tmp_out, expect, shallow=False)
    os.unlink(tmp_out)


def test_load_distance_matrix_file():
    infile = os.path.join(data_dir, "load_distance_matrix_file.txt")
    expect_names = ["sample1", "sample2", "sample3"]
    expect_distances = {(0, 1): 3, (0, 2): 4, (1, 2): 42}
    got_names, got_distances = distances.load_distance_matrix_file(infile)
    assert got_names == expect_names
    assert got_distances == expect_distances

    bad_infile = os.path.join(
        data_dir, "load_distance_matrix_file.wrong_sample_number.txt"
    )
    with pytest.raises(RuntimeError):
        distances.load_distance_matrix_file(bad_infile)
