import numpy as np
from anndata import AnnData  # type: ignore
from beartype.typing import Tuple
from further_types import NormalFormRNASeqData, PositiveInt, NonNegativeInt, Seed



def generate_bulk_data_from_single_cell_data(
        single_cell_data: NormalFormRNASeqData,
        number_of_samples: PositiveInt,
        scale: NonNegativeInt,
        seed: Seed) -> Tuple[NormalFormRNASeqData, AnnData]:
    """Generate bulk sata from single cell data by randomly uniformly set cell numbers for each cell state and
    multiplying it with the single cell matrix"""
    np.random.seed(seed)
    cell_numbers_of_cell_states = AnnData(np.array(np.round(np.random.exponential(
        scale=scale, size=(single_cell_data.n_vars, number_of_samples))), dtype=np.uint64))
    cell_numbers_of_cell_states.var_names = ["simulated_sample_{}".format(i) for i in range(number_of_samples)]
    cell_numbers_of_cell_states.obs_names = single_cell_data.var_names

    bulk_data = AnnData(np.array(np.matmul(single_cell_data.X, cell_numbers_of_cell_states.X), dtype=np.uint64))
    bulk_data.obs_names = single_cell_data.obs_names
    bulk_data.var_names = cell_numbers_of_cell_states.var_names

    return bulk_data, cell_numbers_of_cell_states
