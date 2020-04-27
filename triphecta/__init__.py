from pkg_resources import get_distribution

try:
    __version__ = get_distribution("triphecta").version
except:
    __version__ = "local"


__all__ = [
    "distances",
    "genotypes",
    "phenotypes",
    "phenotype_compare",
    "sample_neighbours_finding",
    "strain_triple",
    "strain_triples",
    "tasks",
    "tree",
    "triphenotops",
    "utils",
    "variant_counts",
    "vcf",
]

from triphecta import *
