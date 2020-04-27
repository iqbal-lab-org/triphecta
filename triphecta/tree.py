import logging

import dendropy

from triphecta import utils


def newick_from_dist_matrix(infile, outfile, method):
    logging.info(f"Loading distance matrix file {infile}")
    with utils.open_file(infile) as f:
        # triphecta saves distance matrix in the standard "phylip" format.
        # First line the number of samples. There is no line of just sample
        # names. This means we need to skip the first line, and then tell
        # dendropy that the first line is not sample names.
        next(f)
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=f, is_first_row_column_names=False, delimiter="\t"
        )

    if method == "upgma":
        logging.info("Calculating upgma tree")
        tree = pdm.upgma_tree()
    elif method == "nj":
        logging.info("Calculating nj tree")
        tree = pdm.nj_tree()
    else:
        raise ValueError(
            f"Got method {method}, but must be upgma or nj. Cannot continue"
        )

    logging.info(f"Writing tree to file {outfile}")
    with utils.open_file(outfile, "w") as f:
        print(
            tree.as_string("newick", suppress_rooting=True).replace("'", ""),
            end="",
            file=f,
        )
    logging.info(f"Finished making tree")
