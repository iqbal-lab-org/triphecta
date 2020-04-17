import collections
import csv

from triphecta import utils

VariantCounts = collections.namedtuple(
    "VariantCounts", ["het", "hom", "null", "het_to_hom"]
)


def save_variant_count_list_to_tsv(var_list, outfile):
    with utils.open_file(outfile, "w") as f:
        print(*VariantCounts._fields, sep="\t", file=f)
        for v in var_list:
            print(*v, sep="\t", file=f)


def load_variant_count_list_from_tsv(infile):
    variants = []
    with utils.open_file(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for d in reader:
            for key in d:
                d[key] = int(d[key])
            variants.append(VariantCounts(**d))
    return variants
