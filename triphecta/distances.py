import csv
import itertools
import multiprocessing
import pickle

import numpy as np

from triphecta import utils, vcf

global vcf_data

# This ended up here so multiprocessing works. vcf_data is a global variable,
# so that no copies of it are made. It is likely to be huge. This function uses
# it read-only, so is ok to have multiple processes all using it.
def _dist_two_samples2(i, j):
    global vcf_data
    return (
        i,
        j,
        np.sum(
            vcf_data[i][0] * vcf_data[j][0] * (vcf_data[i][0] - vcf_data[j][0]) != 0
        ),
    )


def distances_between_vcf_files(
    filenames, threads=1, only_use_pass=True, numeric_filters=None
):
    # filenames = dict of sample -> VCF file
    global vcf_data
    sample_names = []
    vcf_files = []
    for sample_name, vcf_file in filenames.items():
        sample_names.append(sample_name)
        vcf_files.append(vcf_file)

    vcf_data = vcf.load_vcf_files_for_distance_calc(
        vcf_files,
        threads=threads,
        only_use_pass=only_use_pass,
        numeric_filters=numeric_filters,
    )

    with multiprocessing.Pool(processes=threads) as p:
        distance_list = p.starmap(
            _dist_two_samples2, itertools.combinations(range(len(vcf_files)), 2)
        )

    variant_counts = {sample_names[i]: x[1] for i, x in enumerate(vcf_data)}
    dists = {}

    for i, j, dist in distance_list:
        dists[tuple(sorted([sample_names[i], sample_names[j]]))] = dist

    return dists, variant_counts


def pickle_distances_between_vcf_files(
    file_of_vcf_filenames,
    pickle_out,
    threads=1,
    only_use_pass=True,
    numeric_filters=None,
):
    vcf_files = utils.load_file_of_vcf_filenames(
        file_of_vcf_filenames, check_vcf_files_exist=False
    )
    dists, variant_counts = distances_between_vcf_files(
        vcf_files, threads=1, only_use_pass=True, numeric_filters=None
    )
    with open(pickle_out, "wb") as f:
        pickle.dump((dists, variant_counts), f)


def _load_one_sample_distances_file(filename):
    """Loads a distance file into memory. Returns a list of tuples,
       where each tuple is (sample_name, distance)"""
    expect_cols = {"sample", "distance"}
    distances = []

    with open(filename) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not expect_cols.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Error reading distances file {filename}. Expected column names: {','.join(expect_cols)}. Got column names: {','.join(reader.fieldnames)}"
            )
        for row in reader:
            distances.append((row["sample"], float(row["distance"])))

    return distances


def _update_distances_for_one_sample(sample_name, new_distances, all_distances):
    """Updates all distance data in dictionary all_dsitnaces.
    new_distances=list of tuples, made by _load_one_sample_distances_file"""
    for other_sample, distance in new_distances:
        if other_sample == sample_name:
            continue
        key = tuple(sorted([sample_name, other_sample]))
        if key in all_distances and all_distances[key] != distance:
            raise RuntimeError(
                f"Pair of samples seen twice when loading distances, with different distances: {key}"
            )
        all_distances[key] = distance


def _load_sample_distances_file_of_filenames(infile):
    data = []
    expect_cols = {"sample", "distance_file"}
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not expect_cols.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Error reading distances file of filenames {infile}. Expected column names: {','.join(expect_cols)}. Got column names: {','.join(reader.fieldnames)}"
            )
        for row in reader:
            data.append((row["sample"], row["distance_file"]))

    return data


def load_all_one_sample_distances_files(file_of_filenames, threads=1):
    """Loads data from all per sample distances files.
    filenames = dict of sample name -> distance file name.
    <threads> files in parallel"""
    names_and_dist_files = _load_sample_distances_file_of_filenames(file_of_filenames)
    all_distances = {}
    # To reduce memory, load in a batch of files in parallel, then
    # update the distance data in serial.
    for i in range(0, len(names_and_dist_files), threads):
        filenames = [x[1] for x in names_and_dist_files[i : i + threads]]
        with multiprocessing.Pool(processes=threads) as p:
            dist_data = p.map(_load_one_sample_distances_file, filenames)

        for (sample_name, _), new_distances in zip(
            names_and_dist_files[i : i + threads], dist_data
        ):
            _update_distances_for_one_sample(sample_name, new_distances, all_distances)

    return all_distances


def pickle_load_all_one_sample_distances_files(
    file_of_filenames, pickle_out, threads=1
):
    dists = load_all_one_sample_distances_files(file_of_filenames)
    with open(pickle_out, "wb") as f:
        pickle.dump((dists, {}), f)


def load_from_pickle(pickle_file):
    with open(pickle_file, "rb") as f:
        return pickle.load(f)
