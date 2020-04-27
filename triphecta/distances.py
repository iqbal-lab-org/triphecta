import csv
import itertools
import logging
import multiprocessing

import dendropy
import numpy as np

from triphecta import utils, variant_counts, vcf

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
    vcf_names_tsv,
    outprefix,
    threads=1,
    only_use_pass=True,
    numeric_filters=None,
    het_to_hom_key="COV",
    het_to_hom_min_pc_depth=90.0,
    mask_bed_file=None,
):
    logging.info(f"Loading file of VCF filenames {vcf_names_tsv}")
    filenames = utils.load_file_of_vcf_filenames(
        vcf_names_tsv, check_vcf_files_exist=False
    )
    sample_names, vcf_files = zip(*sorted(filenames.items()))
    sample_names = list(sample_names)
    logging.info(f"Found {len(filenames)} VCF files to load")
    logging.info("Getting genotypes from VCF files")

    global vcf_data
    vcf_data = vcf.load_vcf_files_for_distance_calc(
        vcf_files,
        threads=threads,
        only_use_pass=only_use_pass,
        numeric_filters=numeric_filters,
        het_to_hom_key=het_to_hom_key,
        het_to_hom_min_pc_depth=het_to_hom_min_pc_depth,
        mask_bed_file=mask_bed_file,
    )
    logging.info("Finished loading genotypes. Calculating distance matrix")

    with multiprocessing.Pool(processes=threads) as p:
        distance_list = p.starmap(
            _dist_two_samples2, itertools.combinations(range(len(vcf_files)), 2)
        )

    var_counts = [x[1] for x in vcf_data]
    dists = {}
    for i, j, dist in distance_list:
        dists[tuple(sorted([i, j]))] = dist

    logging.info("Finished calculating distance matrix")
    matrix_file = f"{outprefix}.distance_matrix.txt.gz"
    write_distance_matrix_file(sample_names, dists, matrix_file)
    logging.info(f"Saved distance matrix to file {matrix_file}")
    var_counts_file = f"{outprefix}.variant_counts.tsv.gz"
    variant_counts.save_variant_count_list_to_tsv(var_counts, var_counts_file)
    logging.info(f"Saved variant counts file {var_counts_file}")
    return sample_names, dists, var_counts


def _load_one_sample_distances_file(filename):
    """Loads a distance file into memory. Returns a list of tuples,
       where each tuple is (sample_name, distance)"""
    expect_cols = {"sample", "distance"}
    distances = []

    with utils.open_file(filename) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not expect_cols.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Error reading distances file {filename}. Expected column names: {','.join(expect_cols)}. Got column names: {','.join(reader.fieldnames)}"
            )
        for row in reader:
            distances.append((row["sample"], float(row["distance"])))

    return distances


def _update_distances_for_one_sample(
    sample_index, new_distances, all_distances, sample_name_to_index
):
    """Updates all distance data in dictionary all_distances.
    new_distances=list of tuples, made by _load_one_sample_distances_file"""
    for other_sample, distance in new_distances:
        other_index = sample_name_to_index[other_sample]
        if other_index == sample_index:
            continue
        key = tuple(sorted([sample_index, other_index]))
        if key in all_distances and all_distances[key] != distance:
            raise RuntimeError(
                f"Pair of samples seen twice when loading distances, with different distances: {key}"
            )
        all_distances[key] = distance


def _load_sample_distances_file_of_filenames(infile):
    sample_names = []
    distance_files = []
    expect_cols = {"sample", "distance_file"}
    with utils.open_file(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not expect_cols.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Error reading distances file of filenames {infile}. Expected column names: {','.join(expect_cols)}. Got column names: {','.join(reader.fieldnames)}"
            )
        for row in reader:
            sample_names.append(row["sample"])
            distance_files.append(row["distance_file"])

    return sample_names, distance_files


def distances_from_all_one_sample_distances_files(
    file_of_filenames, outfile, threads=1
):
    """Loads data from all per sample distances files.
    filenames = dict of sample name -> distance file name.
    <threads> files in parallel. Writes ditance matrix to outfile, returns
    tuple: sample names list, distances dictionary"""
    sample_names, distance_files = _load_sample_distances_file_of_filenames(
        file_of_filenames
    )
    sample_name_to_index = {name: i for i, name in enumerate(sample_names)}
    all_distances = {}
    # To reduce memory, load in a batch of files in parallel, then
    # update the distance data in serial.
    for i in range(0, len(sample_names), threads):
        end_index = i + threads
        with multiprocessing.Pool(processes=threads) as p:
            dist_data = p.map(
                _load_one_sample_distances_file, distance_files[i:end_index]
            )

        for sample_index, new_distances in zip(range(i, end_index, 1), dist_data):
            _update_distances_for_one_sample(
                sample_index, new_distances, all_distances, sample_name_to_index
            )

    write_distance_matrix_file(sample_names, all_distances, outfile)
    return sample_names, all_distances


def write_distance_matrix_file(sample_names, distance_matrix, outfile):
    with utils.open_file(outfile, "w") as f:
        print(len(sample_names), file=f)
        for i, sample in enumerate(sample_names):
            out = []

            for j, sample2 in enumerate(sample_names):
                if i == j:
                    out.append(0)
                else:
                    out.append(distance_matrix[tuple(sorted([i, j]))])

            print(sample, *out, sep="\t", file=f)


def load_distance_matrix_file(infile):
    sample_names = []
    distances = {}

    with utils.open_file(infile) as f:
        for line_number, line in enumerate(f):
            if line_number == 0:
                try:
                    number_of_samples = int(line.rstrip())
                except:
                    raise RuntimeError(
                        f"Expected first line of distance matrix to contain a number only. Got this: {line}"
                    )

                sample_names = []
            elif line_number == 1:
                sample_names.append(line.split()[0])
                continue
            else:
                fields = line.rstrip().split("\t", maxsplit=line_number)
                sample_names.append(fields[0])
                for i in range(1, line_number):
                    distances[tuple(sorted([line_number - 1, i - 1]))] = float(
                        fields[i]
                    )

    if len(sample_names) != number_of_samples:
        raise RuntimeError(
            f"Expected {number_of_samples} samples in distance matrix file, but got {len(sample_names)}"
        )

    return sample_names, distances


def newick_from_dist_matrix(infile, outfile, method):
    logging.info(f"Loading distance matrix file {infile}")
    with utils.open_file(infile) as f:
        # triphecta saves distance matrix in the standard "phylip" format.
        # First line the number of samples. There is no line of just sample
        # names. This means we need to skip the first line, and then tell
        # dendropy that the first line is not sample names.
        next(f)
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=f, is_first_row_column_names=False, delimiter="\t"
        )

    if method == "upgma":
        logging.info("Calculating upgma tree")
        tree = pdm.upgma_tree()
    elif method == "nj":
        logging.info("Calculating nj tree")
        tree = pdm.nj_tree()
    else:
        raise ValueError(
            f"Got method {method}, but must be upgma or nj. Cannot continue"
        )

    logging.info(f"Writing tree to file {outfile}")
    with utils.open_file(outfile, "w") as f:
        print(
            tree.as_string("newick", suppress_rooting=True).replace("'", ""),
            end="",
            file=f,
        )
    logging.info(f"Finished making tree")
