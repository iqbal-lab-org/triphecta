from triphecta import distances


def run(options):
    sample_names, dists, _ = distances.load_from_pickle(options.distances_file)
    distances.write_distance_matrix_file(sample_names, dists, options.outfile)
