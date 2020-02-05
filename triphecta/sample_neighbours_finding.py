from collections import namedtuple
from operator import attrgetter


RankData = namedtuple(
    "RankData",
    ["sample", "rank_sum", "geno_rank", "pheno_rank", "geno_dist", "pheno_dist"],
)


def _dict_to_value_ranks(d):
    unique_values = sorted(set(d.values()))
    return {x: i for i, x in enumerate(unique_values)}


def _geno_and_pheno_distances_for_one_sample(
    genos, phenos, pheno_compare, sample, top_n_genos=None
):
    phenotype = phenos[sample]
    geno_distances_to_consider = genos.distance_dict(sample, top_n=top_n_genos)
    geno_distances = {}
    pheno_distances = {}

    for other_sample in genos.sample_names():
        if other_sample == sample or other_sample not in geno_distances_to_consider:
            continue

        other_phenotype = phenos[other_sample]
        if not pheno_compare.satisfy_required_differences(phenotype, other_phenotype):
            continue

        geno_distances[other_sample] = geno_distances_to_consider[other_sample]
        pheno_distances[other_sample] = pheno_compare.differences(
            phenotype, other_phenotype
        )

    return geno_distances, pheno_distances


def _geno_and_pheno_distances_to_rank_table(geno_distances, pheno_distances):
    geno_ranks = _dict_to_value_ranks(geno_distances)
    pheno_ranks = _dict_to_value_ranks(pheno_distances)
    rank_table = []

    for s in geno_distances:
        geno_dist = geno_distances[s]
        pheno_dist = pheno_distances[s]
        geno_rank = geno_ranks[geno_dist]
        pheno_rank = pheno_ranks[pheno_dist]

        rank_table.append(
            RankData(
                sample=s,
                rank_sum=geno_rank + pheno_rank,
                geno_rank=geno_rank,
                pheno_rank=pheno_rank,
                geno_dist=geno_dist,
                pheno_dist=pheno_dist,
            )
        )

    rank_table.sort(key=attrgetter("rank_sum", "geno_rank", "pheno_rank"))
    return rank_table


def ranked_neighbours_for_one_sample(
    genos, phenos, pheno_compare, sample, top_n_genos=None
):
    geno_dist, pheno_dist = _geno_and_pheno_distances_for_one_sample(
        genos, phenos, pheno_compare, sample, top_n_genos=top_n_genos
    )
    return _geno_and_pheno_distances_to_rank_table(geno_dist, pheno_dist)
