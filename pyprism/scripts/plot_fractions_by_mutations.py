from setuptools.command.rotate import rotate

import pyprism
import anndata as ad
import os
import matplotlib.pyplot as plt


def main(project: str, gene_id: str) -> None:
    fractions_mutated = ad.read_h5ad(
        os.path.join("..", "..", "fractions_by_mutations", project, gene_id, "fractions_mutated.h5ad"))
    fractions_nonmutated = ad.read_h5ad(
        os.path.join("..", "..", "fractions_by_mutations", project, gene_id, "fractions_nonmutated.h5ad"))
    reference = ad.read_h5ad(
        os.path.join("..", "..", "references_ensembl", "{}.h5ad".format(project)))

    print("Reference")
    #print(reference.uns["cell_state_to_cell_type_dict"])

    print("Mutated")
    print(fractions_mutated)
    fractions_mutated = pyprism.utils.sum_data_parts(fractions_mutated, reference.uns["cell_state_to_cell_type_dict"])
    mean_mutated = fractions_mutated.X.mean(axis=1)
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches(12, 8, forward=True)
    fig.suptitle("{} - {} - mutated".format(project, gene_id))
    plt.bar(fractions_mutated.obs_names, mean_mutated, color="blue")
    plt.xticks(rotation=90)
    plt.ylim(0.0, 0.5)
    plt.show()

    print("Nonmutated")
    print(fractions_nonmutated)
    fractions_nonmutated = pyprism.utils.sum_data_parts(fractions_nonmutated, reference.uns["cell_state_to_cell_type_dict"])
    mean_nonmutated = fractions_nonmutated.X.mean(axis=1)
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches(12, 8, forward=True)
    fig.suptitle("{} - {} - non-mutated".format(project, gene_id))
    plt.bar(fractions_mutated.obs_names, mean_nonmutated, color="blue")
    plt.xticks(rotation=90)
    plt.ylim(0.0, 0.5)
    plt.show()



"""    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches(7 * ncols, 15, forward=True)
    fig.set_dpi(100)
    fig.suptitle("Deconvolution Results", fontsize=20)
    gs = fig.add_gridspec(2, ncols, hspace=0.2, wspace=0.05)
    return gs.subplots(sharey='row')  # sharex='col', sharey='row'
    axes = setup_deconvolution_fig(ncols=bulk_data.n_vars)
    plot_bar_charts_from_anndata(adata=cell_state_fractions, axes=axes[0, :])
    plot_bar_charts_from_anndata(adata=cell_type_fractions, axes=axes[1, :])
    plt.savefig(os.path.join(plot_dir, "{}-{}.png".format(project_name, timestamp)))
    plt.show()"""


if __name__ == "__main__":
    main(project="LUAD", gene_id="ENSG00000133703")
