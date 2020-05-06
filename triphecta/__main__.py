#!/usr/bin/env python3

import argparse
import logging

import triphecta


def main(args=None):
    triphecta_ascii = r"""  _______   _       _               _
 |__   __| (_)     | |             | |
    | |_ __ _ _ __ | |__   ___  ___| |_ __ _
    | | '__| | '_ \| '_ \ / _ \/ __| __/ _` |
    | | |  | | |_) | | | |  __/ (__| || (_| |
    |_|_|  |_| .__/|_| |_|\___|\___|\__\__,_|
             | |
             |_|
"""
    parser = argparse.ArgumentParser(
        prog="triphecta",
        usage="triphecta <command> <options>",
        description=triphecta_ascii,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--version", action="version", version=triphecta.__version__)
    parser.add_argument("--debug", help="Debug mode", action="store_true")
    parser.add_argument("--triphenotops", action="store_true", help=argparse.SUPPRESS)

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ------------------------ vcfs_to_names ----------------------------------
    subparser_vcfs_to_names = subparsers.add_parser(
        "vcfs_to_names",
        help="Get sample names from VCF files",
        usage="triphecta vcfs_to_names [options] <file_of_vcf_filenames> <out.tsv>",
        description="Gets sample names from VCF files, makes TSV file for input to triphecta 'distance_matrix vcf' and to 'triples'",
    )

    subparser_vcfs_to_names.add_argument(
        "--threads",
        type=int,
        help="Number of VCF files to read in parallel [%(default)s]",
        default=1,
        metavar="INT",
    )

    subparser_vcfs_to_names.add_argument(
        "file_of_vcf_filenames", help="File of VCF filenames, one name per line"
    )

    subparser_vcfs_to_names.add_argument("out_tsv", help="Name of output TSV file")

    subparser_vcfs_to_names.set_defaults(func=triphecta.tasks.vcfs_to_names.run)

    # ----------------------- distance_matrix ---------------------------------
    subparser_distance_matrix = subparsers.add_parser(
        "distance_matrix",
        help="Make distance matrix from VCFs or from pairwise pre-made distance files",
        usage="triphecta distance_matrix [options] <vcf|premade> <filenames_tsv> <out>",
        description="Calculates distance between genomes using VCF files, or loads pre-made distances. Saves distance matrix in phylip format",
    )

    subparser_distance_matrix.add_argument(
        "method",
        choices=["vcf", "premade"],
        help="Method to use. Either calculate from VCF files, or use pre-made distances",
    )

    subparser_distance_matrix.add_argument(
        "filenames_tsv",
        help="Name of input data TSV file. Must have 'sample' column, and either 'vcf_file' or 'distance_file' column, depending on method",
    )

    subparser_distance_matrix.add_argument(
        "out",
        help="If method=vcf, prefix of output files. If method=premade, name of single output file",
    )

    subparser_distance_matrix.add_argument(
        "--threads",
        type=int,
        help="Number of threads to use [%(default)s]",
        default=1,
        metavar="INT",
    )

    subparser_distance_matrix.add_argument(
        "--mask_bed_file",
        help="BED file of regions to ignore from all VCF files (tab-delimited, 3 columns: ref_name start end, coords 0-based and end coord not included)",
        metavar="FILENAME",
    )

    subparser_distance_matrix.add_argument(
        "--vcf_ignore_filter_pass",
        action="store_true",
        help="By default, only use VCF records with PASS in the FILTER column, and count all others as a null genotype. Using this flag ignores the filter column, effectively treating every line of the VCF as if it has PASS",
    )

    subparser_distance_matrix.add_argument(
        "--vcf_numeric_filter",
        "-x",
        action="append",
        dest="vcf_numeric_filter",
        help="Add a filter to the VCF parsing, of the form tag_name:(min|max):cutoff. eg GT_CONF:min:10 would count the genotype as null where GT_CONF<10. This option can be used more than once to add more filters",
        metavar="STRING",
    )

    subparser_distance_matrix.add_argument(
        "--het_to_hom_key",
        help="Using this means trying to convert heterozygous calls to homozygous, instead of ignoring them. Use this option to provide the Key to use for allele depths in VCF. It is expected that the corresponding value should be a comma-separated list of allele depths (eg DP4). See also --het_to_hom_cutoff",
        metavar="STRING",
    )

    subparser_distance_matrix.add_argument(
        "--het_to_hom_cutoff",
        type=float,
        help="Only used if --het_to_hom_key is used. --het_to_hom_cutoff X means that a heterozygous call is converted to homozygous for the first allele A found satisfying X <= (100 * depth of A) / (total depth). Otherwise, call is still counted as heterozygous. [%(default)s]",
        default=90.0,
        metavar="FLOAT",
    )

    subparser_distance_matrix.set_defaults(func=triphecta.tasks.distance_matrix.run)

    # ------------------------------ tree -------------------------------------
    subparser_tree = subparsers.add_parser(
        "tree",
        help="Make a tree from distance matrix",
        usage="triphecta tree [options] <nj|upgma> <distance_matrix> <out>",
        description="Make a tree in newick format from distance matrix in phylip format",
    )

    subparser_tree.add_argument(
        "--dendropy",
        action="store_true",
        help="Use dendropy. By default, quicktree is used if it is found in your PATH, otherwise dendropy is used. This option forces dendropy.",
    )

    subparser_tree.add_argument(
        "method", choices=["nj", "upgma"], help="Method to use, either nj or upgma"
    )

    subparser_tree.add_argument(
        "distance_matrix", help="Distance matrix file in phylip format"
    )

    subparser_tree.add_argument("out", help="Name of output newick file")

    subparser_tree.set_defaults(func=triphecta.tasks.tree.run)

    # --------------------- pheno_constraints_template ------------------------
    subparser_pheno_constraints_template = subparsers.add_parser(
        "pheno_constraints_template",
        help="Make a template phenotype constraints file, for use with 'triphecta triples'",
        usage="triphecta pheno_constraints_template <phenos.tsv> <out.json>",
        description="Make a template phenotype constraints file, for use with 'triphecta triples'",
    )

    subparser_pheno_constraints_template.add_argument(
        "phenos_tsv", help="Name of phenotypes TSV file"
    )

    subparser_pheno_constraints_template.add_argument(
        "json_out", help="Name of output phenotypes constraints template JSON file"
    )

    subparser_pheno_constraints_template.set_defaults(
        func=triphecta.tasks.pheno_constraints_template.run
    )

    # ------------------------ find_cases -------------------------------------
    subparser_find_cases = subparsers.add_parser(
        "find_cases",
        help="Find cases, by looking for samples matching certain phenotypes",
        usage="triphecta find_cases [options] <phenos.tsv> <pheno_constraints_json> <outfile>",
        description="Find cases, by looking for samples matching certain phenotypes",
    )

    subparser_find_cases.set_defaults(func=triphecta.tasks.find_cases.run)

    subparser_find_cases.add_argument(
        "--wanted_pheno",
        "-w",
        help="REQUIRED. Phenotype of interest and the value. eg: 'Drug_x,Resistant'. This option can be used more than once, and must be used at least once.",
        action="append",
        required=True,
        metavar="Drug,value",
    )

    subparser_find_cases.add_argument("phenos_tsv", help="Name of phenotypes TSV file")

    subparser_find_cases.add_argument(
        "pheno_constraints_json",
        help="Name of phenotypes constraints file in JSON format",
    )

    subparser_find_cases.add_argument(
        "outfile", help="Name of output file of sample names",
    )

    # ------------------------ triples ----------------------------------------
    subparser_triples = subparsers.add_parser(
        "triples",
        help="Find strain triples and report their variants etc",
        usage="triphecta triples [options] <case_names_file> <vcfs_tsv> <distance_matrix> <phenos_tsv> <pheno_constraints_json> <out>",
        description="Find strain triples and report their variants etc",
    )

    subparser_triples.add_argument(
        "case_names_file", help="Name of file of case sample names. One name per line",
    )

    subparser_triples.add_argument(
        "vcfs_tsv",
        help="Name of input data TSV file. Must have 'sample' and 'vcf_file' columns",
    )

    subparser_triples.add_argument(
        "distance_matrix",
        help="Name of distance matrix file, made by 'triphecta distance_matrix'",
    )

    subparser_triples.add_argument("phenos_tsv", help="Name of phenotypes TSV file")

    subparser_triples.add_argument(
        "pheno_constraints_json",
        help="Name of phenotypes constraints file in JSON format",
    )

    subparser_triples.add_argument("out", help="Prefix of output files")

    subparser_triples.add_argument(
        "--var_counts_file",
        help="Variant counts file *.variant_counts.tsv.gz made by 'triphecta distance_matrix' if method was 'vcf'",
        metavar="FILENAME",
    )

    subparser_triples.add_argument(
        "--mask_bed_file",
        help="BED file of regions to ignore from all VCF files (tab-delimited, 3 columns: ref_name start end, coords 0-based and end coord not included)",
        metavar="FILENAME",
    )

    subparser_triples.add_argument(
        "--top_n_genos",
        help="When finding triples, only consider closest n samples in terms of genetic distance [%(default)s]",
        type=int,
        metavar="INT",
        default=10,
    )

    subparser_triples.add_argument(
        "--max_pheno_diffs",
        help="When finding triples, only consider samples with at most this many phenotype differences [%(default)s]",
        type=int,
        metavar="INT",
        default=1,
    )

    subparser_triples.set_defaults(func=triphecta.tasks.triples.run)

    args = parser.parse_args()
    if args.triphenotops:
        triphecta.triphenotops.roar()

    logging.basicConfig(
        format=f"[%(asctime)s triphecta %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log = logging.getLogger()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
