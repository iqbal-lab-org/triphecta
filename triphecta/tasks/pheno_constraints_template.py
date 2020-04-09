from triphecta import phenotypes


def run(options):
    phenos = phenotypes.Phenotypes(options.phenos_tsv)
    phenos.write_template_constraints_json(options.json_out)
