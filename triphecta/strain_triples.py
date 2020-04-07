import logging

from triphecta import phenotypes, sample_neighbours_finding, strain_triple


class StrainTriples:
    def __init__(self, genos, phenos, pheno_compare, max_pheno_diffs=1, top_n_genos=20):
        self.genos = genos
        self.phenos = phenos
        self.pheno_compare = pheno_compare
        self.max_pheno_diffs = (
            max_pheno_diffs
        )  # TODO needs implementing when finding triples
        self.top_n_genos = top_n_genos
        self.triples = []

    def find_strain_triples(self, wanted_phenos):
        # The initial use case for this was to get a sample that is resistant to
        # drug X, and find two closest (in terms of genomic distance) neighbours
        # that are sensitive to drug X. But we may be interested in other phenotypes,
        # eg high or low growth. Or more than one phenotype.
        # In the code, use the terminology "case" to mean has the phenotype (such
        # as resistant to drug X), and "control" to mean does not have the phenotype
        # (such as sensitive to drug X).

        # eg wanted_phenos = {"d1": True, "d2": 42}
        self.triples = []
        wanted_phenos = {
            k: phenotypes.data_lookup.get(v, v) for k, v in wanted_phenos.items()
        }
        wanted_pheno_keys = set(wanted_phenos.keys())

        for sample_name in self.phenos.phenos:
            if sample_name in self.genos.excluded_samples:
                continue

            if not self.pheno_compare.phenos_agree_on_features(
                self.phenos[sample_name], wanted_phenos, wanted_pheno_keys
            ):
                continue

            logging.info(f"Looking for control samples for case sample '{sample_name}'")
            neighbours = sample_neighbours_finding.ranked_neighbours_for_one_sample(
                self.genos,
                self.phenos,
                self.pheno_compare,
                sample_name,
                top_n_genos=self.top_n_genos,
            )

            if len(neighbours) < 2:
                logging.info(
                    f"Not enough ({len(neighbours)}) controls found for case sample '{sample_name}'"
                )
                continue

            logging.info(
                f"Found {len(neighbours)} potential controls for case sample '{sample_name}"
            )
            logging.info(f"Case: {sample_name}. Control1: {neighbours[0]}")
            logging.info(f"Case: {sample_name}. Control2: {neighbours[1]}")
            self.triples.append(
                strain_triple.StrainTriple(sample_name, neighbours[0], neighbours[1])
            )


    @classmethod
    def _write_triples_names_file(cls, triples, outfile):
        with open(outfile, "w") as f:
            print("triple_id", "case", "control1", "control2", sep="\t", file=f)
            for i, triple in enumerate(triples):
                print(i+1, triple.case, triple.control1, triple.control2, sep="\t", file=f)


