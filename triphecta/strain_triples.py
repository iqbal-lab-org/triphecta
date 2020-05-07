from itertools import repeat
import logging
import multiprocessing
import os

from triphecta import sample_neighbours_finding, strain_triple, utils, vcf


def _process_one_triple(
    triple, triple_index, total_triples, expect_variants, vcf_files, mask, root_out
):
    logging.info(f"Processing triple {triple_index+1} of {total_triples}")
    case_vcf = vcf_files[triple.case]
    control1_vcf = vcf_files[triple.control1.sample]
    control2_vcf = vcf_files[triple.control2.sample]
    triple.set_variants(expect_variants)
    triple.load_variants_from_vcf_files(case_vcf, control1_vcf, control2_vcf)
    triple.update_variants_of_interest()
    outfile = os.path.join(root_out, f"{triple_index+1}.tsv")
    triple.write_variants_of_interest_file(outfile, vcf_records_to_mask=mask)
    triple.clear_variant_calls()
    logging.info(f"Finished triple {triple_index+1} of {total_triples}")
    return triple


class StrainTriples:
    def __init__(
        self,
        genos,
        phenos,
        pheno_compare,
        max_pheno_diffs=1,
        top_n_genos=20,
        processes=1,
    ):
        self.genos = genos
        self.phenos = phenos
        self.pheno_compare = pheno_compare
        self.max_pheno_diffs = max_pheno_diffs
        self.top_n_genos = top_n_genos
        self.triples = []
        self.processes = processes

    def find_strain_triples(self, case_sample_names):
        # The initial use case for this was to get a sample that is resistant to
        # drug X, and find two closest (in terms of genomic distance) neighbours
        # that are sensitive to drug X. But we may be interested in other phenotypes,
        # eg high or low growth. Or more than one phenotype.
        # In the code, use the terminology "case" to mean has the phenotype (such
        # as resistant to drug X), and "control" to mean does not have the phenotype
        # (such as sensitive to drug X).
        for sample_name in case_sample_names:
            if sample_name in self.genos.excluded_samples:
                logging.info(f"Case '{sample_name}' excluded. Skipping")
                continue

            logging.info(f"Looking for control samples for case sample '{sample_name}'")
            neighbours = sample_neighbours_finding.ranked_neighbours_for_one_sample(
                self.genos,
                self.phenos,
                self.pheno_compare,
                sample_name,
                top_n_genos=self.top_n_genos,
                max_pheno_dist=self.max_pheno_diffs,
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
    def _write_triples_names_file(cls, triples, phenos, outfile):
        pheno_names = sorted(list(phenos.pheno_types.keys()))

        with utils.open_file(outfile, "w") as f:
            print(
                "triple_id",
                "case",
                "control1",
                "geno_dist1",
                "pheno_dist1",
                "control2",
                "geno_dist2",
                "pheno_dist2",
                *[f"case_{x}\tcontrol1_{x}\tcontrol2_{x}" for x in pheno_names],
                sep="\t",
                file=f,
            )
            for i, triple in enumerate(triples):
                print(
                    i + 1,
                    triple.case,
                    triple.control1.sample,
                    triple.control1.geno_dist,
                    triple.control1.pheno_dist,
                    triple.control2.sample,
                    triple.control2.geno_dist,
                    triple.control2.pheno_dist,
                    *[
                        f"{phenos[triple.case][x]}\t{phenos[triple.control1.sample][x]}\t{phenos[triple.control2.sample][x]}"
                        for x in pheno_names
                    ],
                    sep="\t",
                    file=f,
                )

    @classmethod
    def _write_variants_summary_file(cls, triples, outfile, vcf_records_to_mask=None):
        with utils.open_file(outfile, "w") as f:
            print(
                "variant_id",
                "in_mask",
                "chrom",
                "pos",
                "ref",
                "alt",
                "freq",
                *[f"Triple.{i+1}" for i in range(len(triples))],
                sep="\t",
                file=f,
            )
            for variant_index, variant in enumerate(triples[0].variants):
                if (
                    vcf_records_to_mask is not None
                    and variant.CHROM in vcf_records_to_mask
                    and variant.POS in vcf_records_to_mask[variant.CHROM]
                ):
                    in_mask = 1
                else:
                    in_mask = 0
                in_triples = [
                    (1 if variant_index in t.variant_indexes_of_interest else 0)
                    for t in triples
                ]
                freq = round(sum(in_triples) / len(in_triples), 4)
                print(
                    variant_index + 1,
                    in_mask,
                    variant.CHROM,
                    variant.POS + 1,
                    variant.REF,
                    ",".join(variant.ALTS),
                    freq,
                    *in_triples,
                    sep="\t",
                    file=f,
                )

    def run_analysis(self, case_sample_names, outprefix, mask_file=None):
        self.find_strain_triples(case_sample_names)
        if len(self.triples) == 0:
            logging.info("No strain triples found. Stopping")
            return

        # The VCFs are expected to have the same positions. Use the first
        # VCF to load the variants and get the mask positions
        vcf_file = self.genos.vcf_files[self.triples[0].case]
        logging.info(f"Load variant positions from first VCF file {vcf_file}")
        _, expect_variants = vcf.load_variant_calls_from_vcf_file(vcf_file)
        if mask_file is None:
            mask = None
        else:
            logging.info(f"Loading mask from file {mask_file}")
            mask = vcf.vcf_to_variant_positions_to_mask_from_bed_file(
                vcf_file, mask_file
            )

        file_per_triple_dir = outprefix + ".triples"
        os.mkdir(file_per_triple_dir)

        with multiprocessing.Pool(processes=self.processes) as pool:
            self.triples = pool.starmap(
                _process_one_triple,
                zip(
                    self.triples,
                    range(len(self.triples)),
                    repeat(len(self.triples)),
                    repeat(expect_variants),
                    repeat(self.genos.vcf_files),
                    repeat(mask),
                    repeat(file_per_triple_dir),
                ),
            )

        triple_names_file = outprefix + ".triple_ids.tsv"
        logging.info(f"Writing file of triple and sample ids {triple_names_file}")
        StrainTriples._write_triples_names_file(
            self.triples, self.phenos, triple_names_file
        )

        variants_file = outprefix + ".variants.tsv"
        logging.info(f"Writing file of variants {variants_file}")
        StrainTriples._write_variants_summary_file(
            self.triples, variants_file, vcf_records_to_mask=mask
        )

        return {
            "triples_names_file": triple_names_file,
            "variants_file": variants_file,
            "triples_dir": file_per_triple_dir,
        }
