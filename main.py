import pyprism
import scanpy as sc
from pathlib import Path

subsample = 1000

store = pyprism.store.Store(path=Path("../store"))
adata = pyprism.datasets.WuEtAl2021(store=store).get("GSM5354522")
print("Number of samples: {}".format(adata.n_obs))
adata = adata[~adata.obs_names.duplicated()].copy()
print("Number of non-duplicate samples: {}".format(adata.n_obs))
adata = pyprism.utils.sample_adata(adata, subsample)
print("Number of subsamples: {}".format(adata.n_obs))
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
pyprism.clustering.calc_kmeans_centroids(adata, n_clusters=10)
print(adata)
adata.write_h5ad("tmp.h5ad")
