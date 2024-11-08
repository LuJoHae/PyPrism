import numpy as np
from anndata import AnnData  # type: ignore
from beartype.typing import Tuple
from pyprism.deconvolution import PositiveInt, NonNegativeInt, Seed, NormalFormRNASeqData


def generate_bulk_data_from_single_cell_data(
        single_cell_data: NormalFormRNASeqData,
        number_of_samples: PositiveInt,
        scale: NonNegativeInt,
        seed: Seed = 0) -> Tuple[NormalFormRNASeqData, AnnData]:
    """Generate bulk sata from single cell data by randomly uniformly set cell numbers for each cell state and
    multiplying it with the single cell matrix"""
    np.random.seed(seed)
    cell_numbers_of_cell_states = AnnData(np.array(np.round(np.random.exponential(
        scale=scale, size=(single_cell_data.n_obs, number_of_samples))), dtype=np.uint64))
    cell_numbers_of_cell_states.var_names = ["sample-{:08d}".format(i) for i in range(number_of_samples)]
    cell_numbers_of_cell_states.obs_names = single_cell_data.obs_names

    bulk_data = AnnData(np.array(np.matmul(cell_numbers_of_cell_states.X.T, single_cell_data.X), dtype=np.uint64))
    bulk_data.var_names = single_cell_data.var_names
    bulk_data.obs_names = cell_numbers_of_cell_states.var_names

    return bulk_data, cell_numbers_of_cell_states.T


def generate_toy_single_cell_data(cell_number, gene_number, seed=0):
    np.random.seed(seed)
    adata = AnnData(np.array(np.random.exponential(size=(cell_number, gene_number), scale=2).round(), dtype=np.uint64))
    adata.obs_names = [("cell-{:08d}").format(i) for i in range(cell_number)]
    adata.var_names = [("gene-{:08d}").format(i) for i in range(gene_number)]
    return adata
