from beartype.typing import Annotated, TypeVar
from beartype.vale import Is
from beartype import beartype
from dataclasses import dataclass
from numpy import uint64, ndarray
from anndata import AnnData
from numpy import generic
from h5py import File as FileH5


DType = TypeVar("DType", bound=generic)

type PositiveInt = Annotated[int, Is[lambda x: x > 0]]
type NonNegativeInt = Annotated[int, Is[lambda x: x >= 0]]
type NonNegativeFloat = Annotated[float, Is[lambda x: x >= 0]]
type Seed = Annotated[int, Is[lambda x: 0 <= x < 2**32]]

"""Normal form is defined as having unscaled count data in an numpy ndarray and ensembl gene ids as obs_names 
(not yet assured)."""
type NormalFormRNASeqData = Annotated[
    AnnData,
    Is[lambda adata: type(adata.X) is ndarray],
    Is[lambda adata: adata.X.dtype == uint64]
]


@beartype
@dataclass(frozen=True)
class DeconvolutionResult:
    """Class for storing the results of a deconvolution. Is returned by deconvolution and also keeps track of the used
    parameters of the convolution."""
    expression_tensor: ndarray
    gene_names: ndarray
    sample_names: ndarray
    cell_state_names: ndarray
    number_of_iterations: PositiveInt

    def __str__(self):
        return "DeconvolutionResult(Genes: {}, Cell States: {}, Samples: {}, Iterations: {})".format(
            *self.expression_tensor.shape, self.number_of_iterations)

    def get_cell_state_fractions(self):
        cell_state_fraction = self.expression_tensor.sum(axis=0)
        cell_state_fraction = cell_state_fraction / cell_state_fraction.sum(axis=0)
        cell_state_fraction = AnnData(cell_state_fraction)
        cell_state_fraction.obs_names = self.cell_state_names
        cell_state_fraction.var_names = self.sample_names
        return cell_state_fraction

    def write_h5dr(self, filename: str) -> None:
        with FileH5(filename, "w") as file:
            file.create_dataset(name="expression_tensor", data=self.expression_tensor)
            file.create_dataset(name="gene_names", data=self.gene_names)
            file.create_dataset(name="sample_names", data=self.sample_names)
            file.create_dataset(name="cell_state_names", data=self.cell_state_names)
            file.create_dataset(name="number_of_iterations", data=self.number_of_iterations)
