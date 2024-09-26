import os
import sys

import pyensembl
import pandas as pd
import numpy as np
import anndata as ad
import pyprism
import pyensembl


def deconvolute_mutations(project, min_number_of_samples=100, number_of_iterations: int = 100):
    reference = ad.read_h5ad("../../references_ensembl/{}.h5ad".format(project))
    adata = ad.read_h5ad("../../samples_mutations_unique-vars/{}.h5ad".format(project))
    print(reference)
    print(adata)

    print(10*"=")
    adata, reference = pyprism.utils.reduce_to_common_var(adata.copy(), reference.copy())
    print(reference)
    print(adata)
    adata.var["mutations_count"] = adata.layers["mutations"].sum(axis=0)
    frequently_mutated_genes = adata.var[adata.var["mutations_count"].to_numpy() > min_number_of_samples].index
    print(frequently_mutated_genes)
    for mutated_gene in frequently_mutated_genes:
        print("===== {} =====".format(mutated_gene))
        index = np.where(adata.var.index == mutated_gene)[0]
        mask = adata.layers["mutations"][:, index]
        mutated_samples = adata[mask, :].copy()
        print(mutated_samples)
        print(reference)
        dr = pyprism.deconvolution.deconvolution(bulk=mutated_samples[1, :], reference=reference, n=number_of_iterations)
        print(dr)
        print(dr.to_df())
        print(mutated_samples)

        break

        non_mutated_samples = adata[~mask, :].copy()
        print(non_mutated_samples)
        dr = pyprism.deconvolution.deconvolution(bulk=non_mutated_samples[1, :].copy(), reference=reference.copy(), n=number_of_iterations)
        print(dr)
    return 0


def add_mutation_data(project):
    adata = ad.read_h5ad("../../samples/{}.h5ad".format(project))
    adata.obs["gene_id_version"] = [gene_id.split(".")[1] for gene_id in adata.obs_names]
    adata.obs_names = [gene_id.split(".")[0] for gene_id in adata.obs_names]
    adata.var["sample_barcode"] = adata.var_names
    adata.var_names = [barcode[:12] for barcode in adata.var_names]

    bdata = adata.copy()
    bdata.X = np.full((adata.n_obs, adata.n_vars), False)

    mutations = pd.read_hdf("../../mutations.h5pd", mode="r")
    mutations = mutations[mutations["project"] == project]

    alice, bob = 0, 0
    for row in mutations.itertuples():
        gene_id, barcode = row[2], row[17][:12]
        try:
            bdata[gene_id, barcode] = True
        except:
            pass

    adata.layers["mutations"] = bdata.X
    adata.write_h5ad("../../samples_mutations/{}.h5ad".format(project), compression="gzip")
    print(adata)
    print(adata.X[:5, :5])


def correct_references(project):
    adata = ad.read_h5ad("../../references/{}/{}_ref.h5ad".format(project, project)).T
    print(adata)
    print(adata.var)
    adata.X = adata.X*1e8-1
    genome = pyensembl.EnsemblRelease(release=111, species="human")
    adata.var["gene_name"] = adata.var_names
    print(adata.var_names)

    gene_ids = []
    for gene_name in adata.var_names:
        try:
            potential_gene_ids = genome.gene_ids_of_gene_name(gene_name.split(".")[0])
            if len(potential_gene_ids) >= 1:
                print()
                gene_ids.append(potential_gene_ids[0])
            else:
                print("BOBO")
                gene_ids.append(None)
        except:
            gene_ids.append(None)
    print(gene_ids[:10])
    print(sum(x is not None for x in gene_ids))
    adata.var_names = gene_ids
    print(adata)
    adata = adata[:, [name is not None for name in adata.var_names]]
    print(adata)
    adata = adata.copy()
    adata.write_h5ad("../../references_ensembl/{}.h5ad".format(project), compression="gzip")
    return 0


def remove_nonunique_vars(adata):
    unique_var_names, var_names_counts = np.unique(adata.var_names, return_counts=True)
    adata.var_names_make_unique()
    return adata[:, unique_var_names[np.where(var_names_counts == 1)]]


def fix_obs(adata):
    adata.obs["patient_barcode"] = adata.obs_names
    adata.obs_names = adata.obs["sample_barcode"]
    return adata


if __name__ == "__main__":
    # deconvolute_mutations("BRCA")
