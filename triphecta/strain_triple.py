import logging

from triphecta import utils, vcf


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

    def set_variants(self, variants):
        self.variants = variants

    def load_variants_from_vcf_files(self, case_vcf, control1_vcf, control2_vcf):
        logging.info(f"Loading VCF file {case_vcf}")
        self.variant_calls[
            "case"
        ], self.variants = vcf.load_variant_calls_from_vcf_file(
            case_vcf, expected_variants=self.variants
        )
        logging.info(f"Loading VCF file {control1_vcf}")
        self.variant_calls["control1"], _ = vcf.load_variant_calls_from_vcf_file(
            control1_vcf, expected_variants=self.variants
        )
        logging.info(f"Loading VCF file {control2_vcf}")
        self.variant_calls["control2"], _ = vcf.load_variant_calls_from_vcf_file(
            control2_vcf, expected_variants=self.variants
        )

    def clear_variant_calls(self):
        self.variant_calls = {"case": None, "control1": None, "control2": None}

    @classmethod
    def genotypes_are_of_interest(cls, case, control1, control2):
        return (
            case is not None
            and control1 is not None
            and control2 is not None
            and len(control1.intersection(control2)) > 0
            and (
                len(case.intersection(control1))
                == 0
                == len(case.intersection(control2))
            )
        )

    def update_variants_of_interest(self):
        self.variant_indexes_of_interest = set()

        for i, variant in enumerate(self.variants):
            if StrainTriple.genotypes_are_of_interest(
                self.variant_calls["case"][i],
                self.variant_calls["control1"][i],
                self.variant_calls["control2"][i],
            ):
                self.variant_indexes_of_interest.add(i)

    @classmethod
    def genotype_to_string(cls, geno):
        if geno is None:
            return "./."
        elif len(geno) == 1:
            g = str(list(geno)[0])
            return f"{g}/{g}"
        else:
            return "/".join(sorted([str(x) for x in geno]))

    def write_variants_of_interest_file(self, filename, vcf_records_to_mask=None):
        with utils.open_file(filename, "w") as f:
            print(
                "variant_id",
                "in_mask",
                "chrom",
                "pos",
                "ref",
                "alt",
                "case",
                "control1",
                "control2",
                sep="\t",
                file=f,
            )
            for i, variant in enumerate(self.variants):
                if i in self.variant_indexes_of_interest:
                    if (
                        vcf_records_to_mask is not None
                        and variant.CHROM in vcf_records_to_mask
                        and variant.POS in vcf_records_to_mask[variant.CHROM]
                    ):
                        in_mask = 1
                    else:
                        in_mask = 0

                    print(
                        i + 1,
                        in_mask,
                        variant.CHROM,
                        variant.POS + 1,
                        variant.REF,
                        ",".join(variant.ALTS),
                        StrainTriple.genotype_to_string(self.variant_calls["case"][i]),
                        StrainTriple.genotype_to_string(
                            self.variant_calls["control1"][i]
                        ),
                        StrainTriple.genotype_to_string(
                            self.variant_calls["control2"][i]
                        ),
                        sep="\t",
                        file=f,
                    )
