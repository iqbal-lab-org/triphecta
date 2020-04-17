from triphecta import distances, utils


def run(options):
    if options.method == "vcf":
        numeric_filters = utils.command_line_filter_list_to_dict(
            options.vcf_numeric_filter
        )
        distances.distances_between_vcf_files(
            options.filenames_tsv,
            options.out,
            threads=options.threads,
            only_use_pass=not options.vcf_ignore_filter_pass,
            numeric_filters=numeric_filters,
            het_to_hom_key=options.het_to_hom_key,
            het_to_hom_min_pc_depth=options.het_to_hom_cutoff,
            mask_bed_file=options.mask_bed_file,
        )
    else:
        sample_names, dists = distances.distances_from_all_one_sample_distances_files(
            options.filenames_tsv, options.out, threads=options.threads
        )
