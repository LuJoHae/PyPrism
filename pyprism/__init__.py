"""
pyprism.

A deconvolution framework for bulk RNA-seq count data based on BayesPrism, implementing the improved algorithm of
InstaPrism in Python, integrated into the ScVerse, and interoperable with AnnData and HDF5 files.
"""

from .nomenclature import remove_version_suffix, create_gene_dict, change_gene_names_to_gene_ids
from .plotting import plot_bar_charts_from_anndata
from .utils import sum_data_parts, remove_zero_obs, reduce_to_common_obs
from .io import read_h5dr
from .simulation import generate_bulk_data_from_single_cell_data
from .deconvolution import multi_deconvolution, deconvolution

__version__ = "0.0.1"
__author__ = "Lukas Jonathan Haeuser"
