from anndata import AnnData
from numpy import float32, float64, array, repeat, unique, round, ravel
from sklearn.cluster import KMeans
from pandas import DataFrame, Series
from scanpy.preprocessing import neighbors, pca
from scanpy.tools import leiden


def calc_kmeans_centroids(adata: AnnData, n_clusters: int):
    """
    WARNING: Potentially not fully deterministic!!!
    Tried to make it deterministic by going to 64 bits floats as with 32 bits one get random deviations of one or two eps
    (i.e. machine precisions).
    """
    method_name = "KMeans-n={}".format(n_clusters)
    adata.X = adata.X.astype(dtype=float64)
    estimator = KMeans(n_clusters=n_clusters, init="k-means++", random_state=0, n_init=1)
    adata.obs[method_name + "_cluster"] = estimator.fit_predict(adata.X)
    adata.varm[method_name] = DataFrame(array(estimator.cluster_centers_.T, dtype=float32),
                                        columns=["cluster_{}".format(i) for i in range(n_clusters)],
                                        index=adata.var_names)


def calc_leiden_centroids(adata: AnnData, resolution: float = 1.0, n_comps: int = 100, n_neighbors: int = 15) -> str:
    method_name = "Leiden-resolution={}".format(round(resolution, 3))
    pca(adata, n_comps=n_comps)
    neighbors(adata, n_pcs=n_comps, n_neighbors=n_neighbors)
    leiden(adata, resolution=resolution, random_state=0, key_added=method_name, flavor="igraph")
    cluster_statistics = unique(array(adata.obs[method_name], dtype=int))
    adata.varm[method_name] = DataFrame(
        [ravel(array(adata[adata.obs[method_name] == str(cluster_index)].X.mean(axis=0))) for cluster_index in cluster_statistics],
        index=cluster_statistics, columns=adata.var_names
    ).T
    return method_name
