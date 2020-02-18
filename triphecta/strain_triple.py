import logging

from triphecta import phenotypes, sample_neighbours_finding, vcf


class StrainTriple:
    def __init__(self, case, control1, control2):
        self.case = case
        self.control1 = control1
        self.control2 = control2
        self.variant_calls = {"case": None, "control1": None, "control2": None}
        self.variants = None
        self.variant_indexes_of_interest = set()

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def load_variants_from_vcf_files(self, case_vcf, control1_vcf, control2_vcf):
        logging.info(f"Loading VCF file {case_vcf}")
        self.variant_calls[
            "case"
        ], self.variants = vcf.load_variant_calls_from_vcf_file(case_vcf)
        logging.info(f"Loading VCF file {control1_vcf}")
        self.variant_calls["control1"], _ = vcf.load_variant_calls_from_vcf_file(
            control1_vcf, expected_variants=self.variants
        )
        logging.info(f"Loading VCF file {control2_vcf}")
        self.variant_calls["control2"], _ = vcf.load_variant_calls_from_vcf_file(
            control2_vcf, expected_variants=self.variants
        )

    def update_variants_of_interest(self):
        self.variant_indexes_of_interest = set()

        for i, variant in enumerate(self.variants):
            # We are interested in where the two controls have the same
            # genotype call, and the case genotype call is different.
            # We can't say anything about null calls, so skip those.
            if (
                self.variant_calls["case"][i] == "."
                or self.variant_calls["control1"][i] == "."
                or self.variant_calls["control2"][i] == "."
                or self.variant_calls["control1"][i]
                != self.variant_calls["control2"][i]
                or self.variant_calls["case"] == self.variant_calls["control1"]
            ):
                continue

            self.variant_indexes_of_interest.add(i)


def find_strain_triples(
    wanted_phenos, genos, phenos, pheno_compare, max_pheno_diffs=1, top_n_genos=20
):
    # The initial use case for this was to get a sample that is resistant to
    # drug X, and find two closest (in terms of genomic distance) neighbours
    # that are sensitive to drug X. But we may be interested in other phenotypes,
    # eg high or low growth. Or more than one phenotype.
    # In the code, use the terminology "case" to mean has the phenotype (such
    # as resistant to drug X), and "control" to mean does not have the phenotype
    # (such as sensitive to drug X).

    # eg wanted_phenos = {"d1": True, "d2": 42}
    strain_triples = []
    wanted_phenos = {
        k: phenotypes.data_lookup.get(v, v) for k, v in wanted_phenos.items()
    }
    wanted_pheno_keys = set(wanted_phenos.keys())

    for sample_name in phenos.phenos:
        if sample_name in genos.excluded_samples:
            continue

        if not pheno_compare.phenos_agree_on_features(
            phenos[sample_name], wanted_phenos, wanted_pheno_keys
        ):
            continue

        logging.info(f"Looking for control samples for case sample '{sample_name}'")
        neighbours = sample_neighbours_finding.ranked_neighbours_for_one_sample(
            genos, phenos, pheno_compare, sample_name, top_n_genos=top_n_genos
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
        strain_triples.append(StrainTriple(sample_name, neighbours[0], neighbours[1]))

    return strain_triples
