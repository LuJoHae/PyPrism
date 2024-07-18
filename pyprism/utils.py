from numpy import array
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


def remove_zero_obs(adata: AnnData) -> AnnData:
    return adata[adata.X.sum(axis=1) > 0]
