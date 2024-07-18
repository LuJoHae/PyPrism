from pyensembl import EnsemblRelease  # type: ignore
from anndata import AnnData  # type: ignore
from beartype.typing import List


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


def change_gene_names_to_gene_ids(adata: AnnData, ensembl_release: int, species: str) -> AnnData:
    gene_dict = create_gene_dict(gene_names=adata.obs_names, ensembl_release=ensembl_release, species=species)
    adata = adata[list(gene_dict.keys()), ]
    adata.obs_names = [gene_dict[gene_name] for gene_name in adata.obs_names]
    return adata
