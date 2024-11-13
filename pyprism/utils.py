from beartype.typing import List
from numpy import array, random, repeat, float32, float64
from anndata import AnnData  # type: ignore
from pandas.core.common import random_state
from pyensembl import EnsemblRelease
from numpy.random import default_rng as rng
from sklearn.cluster import KMeans
from sklearn_extra.cluster import KMedoids


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
    return adata[rng(seed=0).choice(a=adata.n_obs, size=n, replace=False), :].copy()


def change_gene_names_to_gene_ids(adata: AnnData, ensembl_release: int, species: str) -> AnnData:
    gene_dict = create_gene_dict(gene_names=adata.obs_names, ensembl_release=ensembl_release, species=species)
    adata = adata[list(gene_dict.keys()), ]
    adata.obs_names = [gene_dict[gene_name] for gene_name in adata.obs_names]
    return adata


def create_gene_dict(gene_names: List[str], ensembl_release: int, species: str) -> dict:
    ensembl = EnsemblRelease(release=ensembl_release, species=species)
    assert len(gene_names) == len(set(gene_names)), "Gene names not unique!"
    result = {}
    for gene_name in gene_names:
        try:
            gene_ids = ensembl.gene_ids_of_gene_name(gene_name)
            if len(gene_ids) != 1:
                continue
            result[gene_name] = gene_ids[0]
        except:
            continue
    assert len(result) == len(set(result)), "Resulting IDs not unique!"
    return result


def remove_version_suffix(gene_ids: List[str]) -> List[str]:
    return [gene_id.split(".")[0] for gene_id in gene_ids]


MITOCHONDRIA_GENES = [
    'MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND4', 'MT-ND5',
    'MT-ND6', 'MT-CYB'
]



