from anndata import AnnData, concat  # type: ignore
from beartype.typing import Any
from numpy import array, arange
from scanpy.preprocessing import pca, neighbors
from scanpy.tools import umap
from matplotlib.pyplot import scatter, show


def plot_bar_charts_from_anndata(adata: AnnData, axes: Any) -> None:
    for i in range(adata.n_vars):
        data = array([adata.X[:, i]] + [layer[:, i] for layer in adata.layers.values()])
        x = arange(data.shape[1])
        width = 0.3  # the width of the bars
        multiplier = 0
        for measurement in data:
            offset = width * multiplier
            _ = axes[i].bar(x + offset, measurement, width)
            axes[i].set_title(adata.var_names[i])
            axes[i].set_xticks(arange(adata.n_obs),
                               labels=adata.obs_names, rotation=90)
            axes[i].yaxis.grid(True)
            multiplier += 1


def plot_umap_with_varm(adata, varm_key, obs_key=None, size1=20, size2=50):
    varm = AnnData(adata.varm[varm_key]).T
    number_of_clusters = varm.n_obs
    adata = concat([adata, varm], axis=0, join="outer")
    pca(adata)
    neighbors(adata)
    umap(adata)
    for cluster_index in range(number_of_clusters):
        scatter(*adata.obsm["X_umap"][adata.obs[varm_key] == str(cluster_index), :].T, s=size1, marker=".")
    scatter(*adata.obsm["X_umap"][-number_of_clusters:, :].T, s=size2, marker="*", color="black")
    show()
