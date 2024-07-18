from anndata import AnnData  # type: ignore
from beartype.typing import Any
from numpy import array, arange


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
