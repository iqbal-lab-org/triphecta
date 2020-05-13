import json
import logging

from triphecta import genotypes, phenotype_compare, phenotypes, strain_triples, utils


def run(options):
    with open(options.case_names_file) as f:
        case_sample_names = [x.rstrip() for x in f]

    genos = genotypes.Genotypes(
        file_of_vcf_filenames=options.vcfs_tsv,
        distance_matrix_file=options.distance_matrix,
        variant_counts_file=options.var_counts_file,
    )

    phenos = phenotypes.Phenotypes(options.phenos_tsv)

    with open(options.pheno_constraints_json) as f:
        pheno_constraints = json.load(f)

    pheno_compare = phenotype_compare.PhenotypeCompare(pheno_constraints)

    triples = strain_triples.StrainTriples(
        genos,
        phenos,
        pheno_compare,
        top_n_genos=options.top_n_genos,
        max_pheno_diffs=options.max_pheno_diffs,
        processes=options.processes,
    )
    triples.run_analysis(
        case_sample_names, options.out, mask_file=options.mask_bed_file
    )
