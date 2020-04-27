from triphecta import tree


def run(options):
    tree.newick_from_dist_matrix(
        options.distance_matrix,
        options.out,
        options.method,
        force_dendropy=options.dendropy,
    )
