"""
pyprism.

A deconvolution framework for bulk RNA-seq count data based on BayesPrism, implementing the improved algorithm of
InstaPrism in Python, integrated into the ScVerse, and interoperable with AnnData and HDF5 files.
"""

# from .nomenclature import remove_version_suffix, create_gene_dict, change_gene_names_to_gene_ids
# from .plotting import plot_bar_charts_from_anndata
# from .utils import sum_data_parts, remove_zero_obs, reduce_to_common_obs
# from .io import read_h5dr
# from .simulation import generate_bulk_data_from_single_cell_data, generate_toy_single_cell_data
# from .deconvolution import multi_deconvolution, deconvolution_parallel, deconvolution_parallel
# from .datasets import file_hash_md5, get_dataset, get_dataset_source_details

from . import nomenclature
from . import plotting
from . import utils
from . import io
from . import simulation
from . import deconvolution
from . import datasets

__version__ = "0.0.2"
__author__ = "Lukas Jonathan Haeuser"
