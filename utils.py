from numpy import array
from anndata import AnnData


def sum_data_parts(data: AnnData, partition_map: dict) -> AnnData:
    parts = partition_map.keys()
    result = AnnData(array([data[partition_map[part]].X.sum(axis=0) for part in parts]))
    result.var_names = data.var_names
    result.obs_names = parts
    return result