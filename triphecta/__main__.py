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
        description="Gets sample names from VCF files, makes TSV file for input to triphecta distances vcf",
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

    # ------------------------ distances --------------------------------------
    subparser_distances = subparsers.add_parser(
        "distances",
        help="Process genomic distance information",
        usage="triphecta distances <vcf|premade> <data_tsv> <outfile>",
        description="Calculates distance between genomes using VCF files, or loads pre-made distances. Saves data as binary file for use with other triphecta tasks",
    )

    subparser_distances.add_argument(
        "method",
        choices=["vcf", "premade"],
        help="Method to use. Either calculate from VCF files, or use pre-made distances",
    )

    subparser_distances.add_argument(
        "filenames_tsv",
        help="Name of input data TSV file. Must have 'sample' column, and either 'vcf_file' or 'distance_file' column, depending on method",
    )

    subparser_distances.add_argument("outfile", help="Name of output file")

    subparser_distances.add_argument(
        "--threads",
        type=int,
        help="Number of threads to use [%(default)s]",
        default=1,
        metavar="INT",
    )

    subparser_distances.add_argument(
        "--mask_bed_file",
        help="BED file of regions to ignore from all VCF files (tab-delimited, 3 columns: ref_name start end, coords 0-based and end coord not included)",
        metavar="FILENAME",
    )

    subparser_distances.add_argument(
        "--vcf_ignore_filter_pass",
        action="store_true",
        help="By default, only use VCF records with PASS in the FILTER column, and count all others as a null genotype. Using this flag ignores the filter column, effectively treating every line of the VCF as if it has PASS",
    )

    subparser_distances.add_argument(
        "--vcf_numeric_filter",
        "-x",
        action="append",
        dest="vcf_numeric_filter",
        help="Add a filter to the VCF parsing, of the form tag_name:(min|max):cutoff. eg GT_CONF:min:10 would count the genotype as null where GT_CONF<10. This option can be used more than once to add more filters",
        metavar="STRING",
    )

    subparser_distances.add_argument(
        "--het_to_hom_key",
        help="Using this means trying to convert heterozygous calls to homozygous, instead of ignoring them. Use this option to provide the Key to use for allele depths in VCF. It is expected that the corresponding value should be a comma-separated list of allele depths (eg DP4). See also --het_to_hom_cutoff",
        metavar="STRING",
    )

    subparser_distances.add_argument(
        "--het_to_hom_cutoff",
        type=float,
        help="Only used if --het_to_hom_key is used. --het_to_hom_cutoff X means that a heterozygous call is converted to homozygous for the first allele A found satisfying X <= (100 * depth of A) / (total depth). Otherwise, call is still counted as heterozygous. [%(default)s]",
        default=90.0,
        metavar="FLOAT",
    )

    subparser_distances.set_defaults(func=triphecta.tasks.distances.run)

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
