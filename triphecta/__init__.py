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
    "tasks",
    "triphenotops",
    "utils",
    "vcf",
]

from triphecta import *
