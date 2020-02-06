import collections
import functools
import logging
import multiprocessing

import numpy as np

Variant = collections.namedtuple("Variant", ["CHROM", "POS", "REF", "ALTS"])
VariantCounts = collections.namedtuple("VariantCounts", ["het", "hom", "null"])


def vcf_line_to_variant_and_gt(line):
    try:
        chrom, pos, _, ref, alt, _, filter_str, info, format_keys, format_values = line.rstrip().split(
            "\t"
        )
    except:
        raise RuntimeError(f"Wrong number of columns in VCF file at this line:\n{line}")

    try:
        variant = Variant(CHROM=chrom, POS=int(pos) - 1, REF=ref, ALTS=alt.split(","))
    except:
        raise RuntimeError("Error parsing the following line of VCF file:\n{line}")

    if filter_str == "PASS":
        if not format_keys.startswith("GT"):
            raise RuntimeError(
                f"Need GT to be first key in FORMAT column at line:\n{line}"
            )
        gt = set(format_values.split(":")[0].split("/"))
        if len(gt) > 1:
            gt = "."
        else:
            gt = int(gt.pop())
    else:
        gt = "."

    return gt, variant


def load_variant_calls_from_vcf_file(infile, expected_variants=None):
    with open(infile) as f:
        sample_name = None
        calls = []
        checking_variants = True

        if expected_variants is None:
            expected_variants = []
            checking_variants = False

        for line in f:
            if line.startswith("##CHROM"):
                sample_name = line.rstrip().split("\t")[-1]
            elif not line.startswith("#"):
                gt, variant = vcf_line_to_variant_and_gt(line)
                calls.append(gt)

                if checking_variants:
                    if len(calls) - 1 >= len(expected_variants):
                        raise RuntimeError(
                            f"Too many variants in VCF file {infile}. Expected {len(expected_variants)} but got at least one more than that, so stopping"
                        )
                    if expected_variants[len(calls) - 1] != variant:
                        raise RuntimeError(
                            f"Mismatch in variant calls. Expected to get {expected_variants[len(calls)]} but got {variant} in file {infile}. Cannot continue"
                        )
                else:
                    expected_variants.append(variant)

        if sample_name is not None:
            raise RuntimeError(
                f"Did not find sample name in VCF file {infile}. Cannot continue"
            )

        if len(expected_variants) != len(calls):
            raise RuntimeError(
                f"Expected {len(expected_variants)} calls in VCF file {infile} but got {len(calls)}"
            )

    return calls, expected_variants


def load_vcf_file_for_distance_calc(infile, only_use_pass=True, numeric_filters=None):
    """Loads VCF file, returning a numpy array of genotypes, of type uint16.
    0 means unknown genotype. >0 means the allele number (where 1=ref, 2=first alt,
    etc).
    Format of numeric_filters is {"key": (bool, N)}.
    eg "GT_CONF": (True, 10) would require a minimum GT_CONF of 10 to use the
    called genotype. Otherwise the genotype is zero"""
    if numeric_filters is None:
        numeric_filters = {}

    data = []
    variant_counts = {"hom": 0, "het": 0, "null": 0}

    with open(infile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if only_use_pass and fields[6] != "PASS":
                variant_counts["null"] += 1
                data.append(0)
                continue

            try:
                info = dict(zip(fields[8].split(":"), fields[9].split(":")))
                genos = set(info["GT"].split("/"))
            except:
                raise RuntimeError(
                    f"Error parsing final two columns of VCF file {infile} at this line:\n{line}"
                )

            fail_filter = False

            for key, filt in numeric_filters.items():
                if key in info:
                    val = float(info[key])
                    if (filt[0] and val < filt[1]) or (not filt[0] and val > filt[1]):
                        fail_filter = True
                        break

            if fail_filter or len(genos) > 1 or "." in genos:
                data.append(0)
                if len(genos) > 1:
                    variant_counts["het"] += 1
                else:
                    variant_counts["null"] += 1
            else:
                variant_counts["hom"] += 1
                data.append(int(genos.pop()) + 1)

    logging.debug(f"loaded {infile}")
    variant_counts = VariantCounts(**variant_counts)
    return np.array(data, dtype=np.uint16), variant_counts


def load_vcf_files_for_distance_calc(
    filenames, threads=1, only_use_pass=True, numeric_filters=None
):
    if numeric_filters is None:
        numeric_filters = {}
    with multiprocessing.Pool(processes=threads) as p:
        return p.map(
            functools.partial(
                load_vcf_file_for_distance_calc,
                only_use_pass=only_use_pass,
                numeric_filters=numeric_filters,
            ),
            filenames,
        )


def sample_name_from_vcf(infile):
    """Gets sample name from VCF (in its #CHROM... line).
    Assumes the VCF file only conatins one sample"""
    logging.debug(f"Getting sample name from VCF file {infile}")
    with open(infile) as f:
        for line in f:
            if line.startswith("#CHROM"):
                name = line.rstrip().split("\t")[-1]
                logging.debug(f"Found sample name '{name}' from VCF file {infile}")
                return name

    raise RuntimeError(f"#CHROM line not found in file {infile}")


def sample_names_tsv_from_vcf_file_of_filenames(infile, outfile, threads=1):
    """Input is a file of VCF file names, one name per line.
    Writes a TSV file with columns sample_name, vcf_file"""
    with open(infile) as f:
        vcf_files = [x.rstrip() for x in f.readlines()]

    logging.debug(
        f"Getting sample names from {len(vcf_files)} VCF files using {threads} thread(s)"
    )
    with multiprocessing.Pool(processes=threads) as p:
        sample_names = p.map(sample_name_from_vcf, vcf_files)

    assert len(vcf_files) == len(sample_names)
    logging.debug(f"Writing sample/vcf TSV file {outfile}")
    with open(outfile, "w") as f:
        print("sample", "vcf_file", sep="\t", file=f)
        for sample, vcf_file in zip(sample_names, vcf_files):
            print(sample, vcf_file, sep="\t", file=f)
