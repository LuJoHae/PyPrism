import os.path

import pandas as pd

import pyprism
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt


def s(project, min_number_of_samples=100, number_of_iterations: int = 100):
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


def fractions_by_most_frequent_mutations(project: str, min_number_of_samples: int = 100) -> np.array:
    adata = ad.read_h5ad("../../samples_mutations_unique-vars/{}.h5ad".format(project))
    adata.var["mutations_count"] = adata.layers["mutations"].sum(axis=0)
    frequently_mutated_genes = adata.var[adata.var["mutations_count"].to_numpy() > min_number_of_samples].index
    print(adata)
    print(frequently_mutated_genes)
    for mutated_gene in frequently_mutated_genes:
        print("===== {} =====".format(mutated_gene))
        index = np.where(adata.var.index == mutated_gene)[0]
        assert len(index) == 1, "Index length != 1"
        index = index[0]
        mask = adata.layers["mutations"][:, index]
        sample_barcode_with_mutation = adata.obs_names[mask]
        sample_barcode_without_mutation = adata.obs_names[~mask]
        assert len(sample_barcode_with_mutation) + len(sample_barcode_without_mutation) == adata.n_obs
        assert set.union(set(sample_barcode_without_mutation), set(sample_barcode_with_mutation)) == set(adata.obs_names)
        assert set.intersection(set(sample_barcode_with_mutation), set(sample_barcode_without_mutation)) == set()

        if not os.path.exists("../../fractions_by_mutations/{}/{}".format(project, mutated_gene)):
            os.mkdir("../../fractions_by_mutations/{}/{}".format(project, mutated_gene))

        fractions = []
        for sample_barcode in sample_barcode_with_mutation:
            dr = ad.read_h5ad("../../deconvolution_results/{}/{}.h5ad".format(project, sample_barcode))
            fractions.append(pyprism.deconvolution.calculate_fractions(dr).rename(columns={"fraction": sample_barcode}))
        fractions = ad.AnnData(pd.concat(fractions, axis=1))
        fractions.write_h5ad("../../fractions_by_mutations/{}/{}/fractions_mutated.h5ad".format(project, mutated_gene))
        print(fractions)

        fractions = []
        for sample_barcode in sample_barcode_without_mutation:
            dr = ad.read_h5ad("../../deconvolution_results/{}/{}.h5ad".format(project, sample_barcode))
            fractions.append(pyprism.deconvolution.calculate_fractions(dr).rename(columns={"fraction": sample_barcode}))
        fractions = ad.AnnData(pd.concat(fractions, axis=1))
        fractions.write_h5ad("../../fractions_by_mutations/{}/{}/fractions_nonmutated.h5ad".format(project, mutated_gene))
        print(fractions)


if __name__ == "__main__":
    fractions_by_most_frequent_mutations("BRCA")




















