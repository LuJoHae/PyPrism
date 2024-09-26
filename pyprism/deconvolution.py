import numpy as np
import logging
from beartype import beartype
from beartype.typing import Any, TypeVar
from numba import njit
from anndata import AnnData
from pandas import DataFrame

from pyprism.types_and_classes import PositiveInt, NormalFormRNASeqData, DeconvolutionResult


# @beartype
# def multi_deconvolution(bulk_data: NormalFormRNASeqData,
#                         single_cell_reference: NormalFormRNASeqData,
#                         number_of_iterations: PositiveInt
#                         ) -> DeconvolutionResult:
#
#     # Normalize columns (cell types)
#     normalized_reference = single_cell_reference.X / single_cell_reference.X.sum(axis=0)
#
#     # Initial guess (replace with  initial construction of expression_tensor)
#     cell_state_fraction = np.ones((single_cell_reference.n_vars, bulk_data.n_vars)) / single_cell_reference.n_vars
#
#     # Iteration scheme
#     for i in range(number_of_iterations):
#         expression_tensor = normalized_reference[:, :, None] * cell_state_fraction[None, :, :]
#         expression_tensor = expression_tensor / expression_tensor.sum(axis=1)[:, None, :]
#         expression_tensor = bulk_data.X[:, None, :] * expression_tensor
#         cell_state_fraction = expression_tensor.sum(axis=0)
#         cell_state_fraction = cell_state_fraction / cell_state_fraction.sum(axis=0)
#
#
#     return DeconvolutionResult(expression_tensor=expression_tensor,
#                                gene_names=bulk_data.obs_names.to_numpy(),
#                                sample_names=bulk_data.var_names.to_numpy(),
#                                cell_state_names=single_cell_reference.var_names.to_numpy(),
#                                number_of_iterations=number_of_iterations)
#
#
# @njit(parallel=True)
# def multi_parallel_deconvolution(bulk_data,
#                                  single_cell_reference,
#                                  number_of_iterations
#                                  ):
#
#     # Normalize columns (cell types)
#     normalized_reference = single_cell_reference / single_cell_reference.sum(axis=0)
#
#     # Initial guess (replace with  initial construction of expression_tensor)
#     cell_state_fraction = np.ones((single_cell_reference.shape[1], bulk_data.shape[1])) / single_cell_reference.shape[1]
#
#     # Iteration scheme
#     for i in range(number_of_iterations):
#         expression_tensor = normalized_reference.reshape((normalized_reference.shape[0], normalized_reference.shape[1], 1)) * cell_state_fraction.reshape((1, cell_state_fraction.shape[0], cell_state_fraction.shape[1]))
#         expression_tensor = expression_tensor / expression_tensor.sum(axis=1).reshape((expression_tensor.shape[0], 1, expression_tensor.shape[2]))
#         expression_tensor = bulk_data.reshape((bulk_data.shape[0], 1, bulk_data.shape[1])) * expression_tensor
#         cell_state_fraction = expression_tensor.sum(axis=0)
#         cell_state_fraction = cell_state_fraction / cell_state_fraction.sum(axis=0)
#
#
#     return expression_tensor
#
#
# @beartype
# def _deconvolution(bulk_data: NormalFormRNASeqData,
#                   single_cell_reference: NormalFormRNASeqData,
#                   number_of_iterations: PositiveInt
#                   ) -> Any:
#
#     # Normalize columns (cell types)
#     normalized_reference = single_cell_reference.X / single_cell_reference.X.sum(axis=0)
#
#     # Initial guess
#     cell_state_fraction = np.ones(single_cell_reference.n_vars) / single_cell_reference.n_vars
#
#     # Iteration scheme
#     for i in range(number_of_iterations):
#         B = normalized_reference * cell_state_fraction[None, :]
#         B = B / B.sum(axis=1)[:, None]
#         B = B * bulk_data.X
#         cell_state_fraction = B.sum(axis=0)
#         cell_state_fraction = cell_state_fraction / cell_state_fraction.sum()
#
#     return cell_state_fraction

@njit(parallel=True)
def _deconvolution(bulk: np.array, reference: np.array, n: np.array, eps: float):
    c, g = reference.shape

    # Initialize expression matrix
    B = np.ones(shape=(c, g)) / np.repeat(c * g, c * g).reshape((c, g))

    # Iteration scheme
    for i in range(n):
        B = reference * np.repeat(B.sum(axis=1), g).reshape((c, g))
        B = B / (np.repeat(B.sum(axis=0), c).reshape((g, c)).T + eps)
        B = B * np.repeat(bulk, c).reshape((g, c)).T

    return B


def single_deconvolution(bulk: AnnData, reference: AnnData, n: int = 10, library_normalization: bool = True):
    obs = reference.obs
    var = reference.var
    reference = np.array(reference.X, dtype=np.float64)
    bulk = np.array(bulk.X, dtype=np.float64)

    # library normalize reference
    if library_normalization:
        c, g = reference.shape
        reference = reference / np.repeat(reference.sum(axis=1), g).reshape((c, g))

    result = _deconvolution(bulk=bulk, reference=reference, n=n, eps=np.finfo(float).eps)
    result = AnnData(result)
    result.obs = obs
    result.var = var

    return result


def calculate_fractions(B: AnnData):
    theta = B.X.sum(axis=1)
    theta = theta / theta.sum()
    theta = DataFrame(theta, index=B.obs_names, columns=["fraction"])
    return theta


def get_cell_type_fractions(reference, cell_type_label, fractions_label):
    cell_types = np.unique(reference.obs[cell_type_label])
    result = {type: 0.0 for type in cell_types}
    for cell_type in cell_types:
        df = reference.obs[reference.obs[cell_type_label] == cell_type]
        result[cell_type] = df[fractions_label].sum()
    return result
