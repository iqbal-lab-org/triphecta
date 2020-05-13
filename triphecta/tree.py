import logging
import shutil

import dendropy

from triphecta import utils


def dendropy_newick_from_dist_matrix(infile, outfile, method):
    logging.info("Calculating tree using dendropy")
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


def quicktree_newick_from_dist_matrix(infile, outfile, method):
    # quicktree can't read gzip file
    if infile.endswith(".gz"):
        infile = f"<(zcat {infile})"

    method_opt = " -upgma " if method == "upgma" else " "
    command = f"quicktree{method_opt}-in m -out t {infile} > {outfile}"
    logging.info("Calculating tree using quicktree.")
    utils.syscall(command)


def newick_from_dist_matrix(infile, outfile, method, force_dendropy=False):
    if force_dendropy or shutil.which("quicktree") is None:
        dendropy_newick_from_dist_matrix(infile, outfile, method)
    else:
        quicktree_newick_from_dist_matrix(infile, outfile, method)

    logging.info("Finished making tree")
