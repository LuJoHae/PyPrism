from anndata import AnnData
from numpy import float32, float64, array, repeat
from sklearn.cluster import KMeans
from pandas import DataFrame


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
    adata.varm[method_name] = DataFrame(array(estimator.cluster_centers_, dtype=float32), columns=adata.var_names,
                                        index=["cluster_{}".format(i) for i in range(n_clusters)])
