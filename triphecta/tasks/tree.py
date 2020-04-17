from triphecta import distances


def run(options):
    distances.newick_from_dist_matrix(
        options.distance_matrix, options.out, options.method
    )
