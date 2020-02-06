from triphecta import vcf


def run(options):
    vcf.sample_names_tsv_from_vcf_file_of_filenames(
        options.file_of_vcf_filenames, options.out_tsv, threads=options.threads
    )
