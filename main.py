import pyprism
import logging

logging.getLogger("pyprism").setLevel(logging.DEBUG)

store = pyprism.store.Store("./store")

data = pyprism.datasets.deconvolution_fractions_tcga_all(store).get()