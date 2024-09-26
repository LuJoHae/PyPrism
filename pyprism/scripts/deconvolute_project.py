import anndata as ad
import numpy as np
import pyprism
import os

if __name__ == '__main__':
    project = "BRCA"
    output_dir = "../../results/BRCA"
    n=100

    reference = ad.read_h5ad("../../references_ensembl/{}.h5ad".format(project))
    bulk = ad.read_h5ad("../../samples_mutations_unique-vars/{}.h5ad".format(project))
    bulk, reference = pyprism.utils.reduce_to_common_var(bulk, reference)

    for i in range(bulk.n_obs):
        sample_name = bulk.obs_names[i]
        print(i, sample_name)
        dr = pyprism.deconvolution.single_deconvolution(bulk=bulk[i, :], reference=reference, n=n)
        dr.write_h5ad(os.path.join(output_dir, "{}.h5ad".format(sample_name)))

