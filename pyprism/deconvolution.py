import numpy as np
import logging
from beartype import beartype
from beartype.typing import Any, TypeVar

from pyprism.types_and_classes import PositiveInt, NormalFormRNASeqData, DeconvolutionResult


@beartype
def multi_deconvolution(bulk_data: NormalFormRNASeqData,
                        single_cell_reference: NormalFormRNASeqData,
                        number_of_iterations: PositiveInt
                        ) -> DeconvolutionResult:

    # Normalize columns (cell types)
    normalized_reference = single_cell_reference.X / single_cell_reference.X.sum(axis=0)

    # Initial guess (replace with  initial construction of expression_tensor)
    cell_state_fraction = np.ones((single_cell_reference.n_vars, bulk_data.n_vars)) / single_cell_reference.n_vars

    # Iteration scheme
    for i in range(number_of_iterations):
        expression_tensor = normalized_reference[:, :, None] * cell_state_fraction[None, :, :]
        expression_tensor = expression_tensor / expression_tensor.sum(axis=1)[:, None, :]
        expression_tensor = bulk_data.X[:, None, :] * expression_tensor
        cell_state_fraction = expression_tensor.sum(axis=0)
        cell_state_fraction = cell_state_fraction / cell_state_fraction.sum(axis=0)


    return DeconvolutionResult(expression_tensor=expression_tensor,
                               gene_names=bulk_data.obs_names.to_numpy(),
                               sample_names=bulk_data.var_names.to_numpy(),
                               cell_state_names=single_cell_reference.var_names.to_numpy(),
                               number_of_iterations=number_of_iterations)


@beartype
def deconvolution(bulk_data: NormalFormRNASeqData,
                  single_cell_reference: NormalFormRNASeqData,
                  number_of_iterations: PositiveInt
                  ) -> Any:

    # Normalize columns (cell types)
    normalized_reference = single_cell_reference.X / single_cell_reference.X.sum(axis=0)

    # Initial guess
    cell_state_fraction = np.ones(single_cell_reference.n_vars) / single_cell_reference.n_vars

    # Iteration scheme
    for i in range(number_of_iterations):
        B = normalized_reference * cell_state_fraction[None, :]
        B = B / B.sum(axis=1)[:, None]
        B = B * bulk_data.X
        cell_state_fraction = B.sum(axis=0)
        cell_state_fraction = cell_state_fraction / cell_state_fraction.sum()

    return cell_state_fraction