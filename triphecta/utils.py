from contextlib import contextmanager
import csv
import gzip
import logging
import os
import subprocess
import sys

from triphecta import phenotypes


def rm_rf(*paths):
    for path in paths:
        subprocess.check_output(f"rm -rf {path}", shell=True)


@contextmanager
def open_file(filename, mode="r"):
    if filename.endswith(".gz"):
        try:
            f = gzip.open(filename, f"{mode}t", compresslevel=9)
        except:
            raise OSError(
                f"Error opening gzip file '{filename}' in mode '{mode}'. Cannot continue"
            )
    elif filename.endswith(".bcf"):
        raise NotImplementedError("Opening BCF files is not implemented")
    else:
        try:
            f = open(filename, mode)
        except:
            raise OSError(
                f"Error opening file '{filename}' in mode '{mode}'. Cannot continue"
            )

    yield f
    f.close()


def load_file_of_vcf_filenames(filename, check_vcf_files_exist=True):
    data = {}
    with open_file(filename) as f:
        reader = csv.DictReader(f, delimiter="\t")
        expect_cols = {"sample", "vcf_file"}
        if not expect_cols.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Error reading file of VCF filenames  {filename}. Expected column names: {','.join(expect_cols)}. Got column names: {','.join(reader.fieldnames)}"
            )
        for row in reader:
            if row["sample"] in data:
                raise RuntimeError(
                    f"Duplicated sample '{row['sample']}' in file of VCF filnenames {filename}. Cannot continue"
                )
            if check_vcf_files_exist and not os.path.exists(row["vcf_file"]):
                raise FileNotFoundError(f"Did not find VCF file '{row['vcf_file']}'")

            data[row["sample"]] = row["vcf_file"]

    return data


def command_line_filter_list_to_dict(filter_list):
    if filter_list is None:
        return {}

    filters = {}

    for filter_string in filter_list:
        try:
            name, min_or_max, cutoff = filter_string.split(":")
            cutoff = float(cutoff)
        except:
            raise RuntimeError(f"Error parsing filter {filter_string}")

        if min_or_max == "min":
            m = True
        elif min_or_max == "max":
            m = False
        else:
            raise RuntimeError(f"Error parsing filter {filter_string}")

        filters[name] = (m, cutoff)

    return filters


def command_line_wanted_phenos_to_dict(pheno_list):
    if pheno_list is None:
        return {}

    wanted_phenos = {}

    for string in pheno_list:
        pheno, value = string.split(",")
        try:
            wanted_phenos[pheno] = phenotypes.Phenotypes.convert_one_variable_string(
                value
            )
        except TypeError:
            raise TypeError(f"Error parsing command line phenotypes string '{string}'")

    return wanted_phenos


def syscall(command):
    logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        executable="/bin/bash",
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logging.info(f"Return code: {completed_process.returncode}")
    if completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "Output from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise RuntimeError("Error in system call. Cannot continue")

    logging.info(f"stdout:\n{completed_process.stdout.rstrip()}")
    logging.info(f"stderr:\n{completed_process.stderr.rstrip()}")
    return completed_process
