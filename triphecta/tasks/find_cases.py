import json

from triphecta import phenotype_compare, phenotypes, utils


def run(options):
    wanted_phenos = utils.command_line_wanted_phenos_to_dict(options.wanted_pheno)
    phenos = phenotypes.Phenotypes(options.phenos_tsv)
    with open(options.pheno_constraints_json) as f:
        pheno_constraints = json.load(f)
    pheno_compare = phenotype_compare.PhenotypeCompare(pheno_constraints)
    samples = phenos.find_matching_cases(wanted_phenos, pheno_compare)
    if len(samples) == 0:
        raise RuntimeError("No matching samples found")
    with open(options.outfile, "w") as f:
        print(*samples, sep="\n", file=f)
