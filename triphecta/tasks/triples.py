import json
import logging

from triphecta import genotypes, phenotype_compare, phenotypes, strain_triples, utils


def run(options):
    genos = genotypes.Genotypes(
        file_of_vcf_filenames=options.vcfs_tsv,
        distances_pickle_file=options.distances_file,
    )

    phenos = phenotypes.Phenotypes(options.phenos_tsv)
    wanted_phenos = utils.command_line_wanted_phenos_to_dict(options.wanted_pheno)

    with open(options.pheno_constraints_json) as f:
        pheno_constraints = json.load(f)

    for drug in wanted_phenos:
        if pheno_constraints[drug]["must_be_same"]:
            logging.info(
                f"Changing 'must_be_same' from true to false for phenotype '{drug}', because it was given by the option --wanted_pheno"
            )
            pheno_constraints[drug]["must_be_same"] = False

    pheno_compare = phenotype_compare.PhenotypeCompare(pheno_constraints)

    triples = strain_triples.StrainTriples(
        genos,
        phenos,
        pheno_compare,
        top_n_genos=options.top_n_genos,
        max_pheno_diffs=options.max_pheno_diffs,
    )
    triples.run_analysis(wanted_phenos, options.out, mask_file=options.mask_file)
