import matplotlib.pyplot as plt
import os.path
from anndata import read_h5ad
from beartype.typing import Any
from datetime import datetime
from deconvolution import multi_deconvolution
from utils import sum_data_parts
from simulation import generate_bulk_data_from_single_cell_data
from plotting import plot_bar_charts_from_anndata


def setup_deconvolution_fig() -> Any:
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches(20, 15, forward=True)
    fig.set_dpi(100)
    fig.suptitle("Deconvolution Results", fontsize=20)
    gs = fig.add_gridspec(2, bulk_data.n_vars, hspace=0.2, wspace=0.05)
    return gs.subplots(sharey='row')  # sharex='col', sharey='row'


if __name__ == "__main__":

    # OS calls
    timestamp = "{:%Y%m%d%H%M%S}".format(datetime.now())
    plot_dir = os.path.join("plots")
    results_dir = os.path.join("results")
    signatures_dir = os.path.join("signatures")
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    if not os.path.exists(plot_dir):
        raise "Signature directory \"{}\" does not exist!".format(signatures_dir)

    # Load data
    reference = read_h5ad(os.path.join("signatures", "ov_ref.h5ad"))
    assert (reference.X >= 0).all(), "Reference has negative entries!"
    reference.X += 1  # add one pseudo count to each entry

    # Simulate bulk data
    bulk_data, cell_numbers_of_cell_states = generate_bulk_data_from_single_cell_data(
        reference, number_of_samples=3, scale=100000, seed=0)

    # Deconvolution
    deconvolution_result = multi_deconvolution(
        bulk_data=bulk_data, single_cell_reference=reference, number_of_iterations=10)

    # Save results to disk
    deconvolution_result.write_h5dr(filename=os.path.join(results_dir, "test-{}.h5dr".format(timestamp)))

    # Calculate fractions
    cell_state_fractions = deconvolution_result.get_cell_state_fractions()
    cell_type_fractions = sum_data_parts(cell_state_fractions, reference.uns["cell_state_to_cell_type_dict"])

    # Add ground truth
    cell_state_fractions.layers["true"] = cell_numbers_of_cell_states.X / cell_numbers_of_cell_states.X.sum(axis=0)
    cell_numbers_of_cell_types = sum_data_parts(
        cell_numbers_of_cell_states, reference.uns["cell_state_to_cell_type_dict"])
    cell_type_fractions.layers["true"] = cell_numbers_of_cell_types.X / cell_numbers_of_cell_types.X.sum(axis=0)

    # Plotting
    axes = setup_deconvolution_fig()
    plot_bar_charts_from_anndata(adata=cell_state_fractions, axes=axes[0, :])
    plot_bar_charts_from_anndata(adata=cell_type_fractions, axes=axes[1, :])
    plt.savefig(os.path.join(plot_dir, "test-{}.png".format(timestamp)))
    plt.show()
