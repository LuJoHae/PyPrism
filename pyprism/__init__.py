"""
pyprism.

A deconvolution framework for bulk RNA-seq count data based on BayesPrism, implementing the improved algorithm of
InstaPrism in Python, integrated into the ScVerse, and interoperable with AnnData and HDF5 files.
"""

from . import (
    nomenclature,
    io,
    plotting,
    simulation,
    utils,
    deconvolution
)

__version__ = "0.0.1"
__author__ = "Lukas Jonathan Haeuser"
