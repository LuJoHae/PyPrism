"""
pyprism.

A deconvolution framework for bulk RNA-seq count data based on BayesPrism, implementing the improved algorithm of
InstaPrism in Python, integrated into the ScVerse, and interoperable with AnnData and HDF5 files.
"""


from . import plotting
from . import utils
from . import simulation
from . import deconvolution
from . import datasets
from . import store
from . import clustering

__version__ = "0.0.2"
__author__ = "Lukas Jonathan Haeuser"

_STORE_DIRECTORY = "./.store"
