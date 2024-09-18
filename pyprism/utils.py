from numpy import array, random
from anndata import AnnData  # type: ignore


def sum_data_parts(data: AnnData, partition_map: dict) -> AnnData:
    parts = partition_map.keys()
    result = AnnData(array([data[partition_map[part]].X.sum(axis=0) for part in parts]))
    result.var_names = data.var_names
    result.obs_names = parts
    return result


def reduce_to_common_obs(adata, bdata):
    obs_list = sorted(list(set.intersection(set(adata.obs_names), set(bdata.obs_names))))
    return adata[obs_list, :].copy(), bdata[obs_list, :].copy()


def reduce_to_common_var(adata, bdata):
    var_list = sorted(list(set.intersection(set(adata.var_names), set(bdata.var_names))))
    return adata[:, var_list].copy(), bdata[:, var_list].copy()


def remove_zero_obs(adata: AnnData) -> AnnData:
    return adata[adata.X.sum(axis=1) > 0]


def sample_adata(adata, n, seed=0):
    random.seed(seed)
    return adata[random.randint(0, adata.n_obs, n), :].copy()