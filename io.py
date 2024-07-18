from further_types import DeconvolutionResult
from h5py import File as FileH5  # type: ignore
from numpy import array


def read_h5dr(filename: str) -> DeconvolutionResult:
    with FileH5(filename, "r") as file:
        return DeconvolutionResult(
            expression_tensor=file["expression_tensor"][:],
            gene_names=file["gene_names"][:],
            sample_names=file["sample_names"][:],
            cell_state_names=file["cell_state_names"][:],
            number_of_iterations=int(array(file["number_of_iterations"]))
        )