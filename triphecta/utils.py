import csv
import os


def load_file_of_vcf_filenames(filename, check_vcf_files_exist=True):
    data = {}
    with open(filename) as f:
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
