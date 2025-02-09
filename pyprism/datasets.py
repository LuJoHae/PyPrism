import logging
import pandas as pd

from pyprism.deconvolution import calculate_fractions
from pyprism.store import Store, StoreElement, StoreElementContainer, Hash, download_gsm, GSM
from pyprism.tcga import patient

from anndata import AnnData, read_h5ad
from anndata import concat as ad_concat
from collections import OrderedDict
from functools import partial
from io import TextIOWrapper
from numpy import repeat, unique, array, zeros, empty, concatenate, ndarray
from os import listdir, mkdir
from pandas import read_hdf as pd_read_hdf, read_hdf
from pandas import concat as pd_concat
from pandas import read_excel
from pathlib import Path
from pickle import dump, load
from re import match
from shutil import copyfile
from xarray import DataArray, Dataset, DataTree, open_dataset, open_dataarray, open_datatree
from xarray import concat as xr_concat

logger = logging.getLogger('pyprism')


class WuEtAl2021(StoreElementContainer):

    def __init__(self, store: Store):
        self.store = store
        self.gsm_list = [GSM(f"GSM53545{i:02d}") for i in range(13, 39)]
        self.darwin_hashes = [
            Hash("d39e99b726e5e5c8a913baf0129c012f"),  # GSM5354513
            Hash("cc1c8ca5f8f61092f809c3e26ea9316f"),  # GSM5354514
            Hash("bc0debfa94905d725bc2fc06f9ba0002"),  # GSM5354515
            Hash("dffea4cc4503b293114beb633e9c35f9"),  # GSM5354516
            Hash("f99ab8b3dd692d224ae20bd7264c6263"),  # GSM5354517
            Hash("497fad0872f7d63a06b488bdb599daef"),  # GSM5354518
            Hash("7d0db0145865d63f088c54e4a2ad2041"),  # GSM5354519
            Hash("187dc4ad7db38315323aa1f6dfe1ee0e"),  # GSM5354520
            Hash("6e25f3167c62d7e46666a391b36b1134"),  # GSM5354521
            Hash("9cc4259c0780ae1c28daa7def42b9bfd"),  # GSM5354522
            Hash("048cf20f2adfa25246e6f49c8f94dc66"),  # GSM5354523
            Hash("c91ff6df70642c650ae046f9639a7efc"),  # GSM5354524
            Hash("036b84c2a907672a575d6e75f5bc78a7"),  # GSM5354525
            Hash("abbef152dca951dc26e9ab4567a6fd5a"),  # GSM5354526
            Hash("366cc6d1c0534dc25e7c6968400ff8a1"),  # GSM5354527
            Hash("409f38ef3cc71e9fa7bb94ddedc3f4b1"),  # GSM5354528
            Hash("db27525e9d2377e57380793919be8ef1"),  # GSM5354529
            Hash("1be71dbd7d795773a5cd3c865d92ec23"),  # GSM5354530
            Hash("7ac4359e4b9b5283cd0c045f6d167e80"),  # GSM5354531
            Hash("3960d720b64a14d30911352ea399f202"),  # GSM5354532
            Hash("621c4804f045888a0da045c336b1d430"),  # GSM5354533
            Hash("22cfd27bd71c151a7e314cd13929dc38"),  # GSM5354534
            Hash("fcaeb44787a792bfe6532f231f7d948d"),  # GSM5354535
            Hash("cd8d36c811e8c90dead0692f71715f4a"),  # GSM5354536
            Hash("a1b383932c1c8cc60861b3ce8e3dd615"),  # GSM5354537
            Hash("8e4dce62cd14ce526bbae8a31656fbfe")   # GSM5354538
        ]
        self.linux_hashes = [
            Hash("bca6fab1b19a441f52d392383ee15c12"),  # GSM5354513
            Hash("0902b8a3e71e7d3946e11337e9092db6"),  # GSM5354514
            Hash("f1a031c2eafa0cce8969d11b87daea06"),  # GSM5354515
            Hash("3d28a4fad7e920758341e2d8c5752418"),  # GSM5354516
            Hash("d2dfa38e45067697dee92ee959cfeac6"),  # GSM5354517
            Hash("eba547ab308198b28f1b4dc722be7fee"),  # GSM5354518
            Hash("079abcc38163434639a3854bcab3a05a"),  # GSM5354519
            Hash("6473e1db9ce0f769fc64241dd0a54168"),  # GSM5354520
            Hash("2d10d93de16f9cd77d099074d76cbd70"),  # GSM5354521
            Hash("afca7b0cbc4ffb5b166c35c2778cf419"),  # GSM5354522
            Hash("a6ef00fc4c661be694d9240343311dbf"),  # GSM5354523
            Hash("1a44090f93577dd2e236f02fbed98560"),  # GSM5354524
            Hash("fedd1dde1f6ec4cd468b218311325aa7"),  # GSM5354525
            Hash("dfb991ee8d050b5b9d41c7e86833223f"),  # GSM5354526
            Hash("d6580ade6fb5ddef67ed7213f4d0cba8"),  # GSM5354527
            Hash("8204aa5649fb8d01aad02b853597f7fd"),  # GSM5354528
            Hash("1cf9f2cf384eec009e50e1c633bc0455"),  # GSM5354529
            Hash("3d2c2bf8af127ba3aa6d9a9384d82b8b"),  # GSM5354530
            Hash("b7856f0c5dd3381f9be638c7ffd2ce93"),  # GSM5354531
            Hash("0edcdbc53c2191a74614eacb11e21062"),  # GSM5354532
            Hash("a286a242956dcdacebd0e02fd8a789d1"),  # GSM5354533
            Hash("1f1fdc3e9e2a00b782bf8ab107c9d632"),  # GSM5354534
            Hash("c0dfb2bced3a939e92efd810a8370031"),  # GSM5354535
            Hash("46f0aca3609bd7f441d893fa3c88525d"),  # GSM5354536
            Hash("80205dc798d15c6c12b01dd03dd4d0b3"),  # GSM5354537
            Hash("38fc4d64c42e9d67f1be44f02ecdd28e")   # GSM5354538
        ]

        def derivation(gsm: GSM, output_dir: TextIOWrapper, derivation_store_files) -> Path:  # TODO: from derivation_store_files to derivation_store_file_HANDLES
            gsm_filepath = download_gsm(gsm=gsm, output_dir=output_dir)
            logging.debug(f"GSM filepath: {gsm_filepath}")
            return gsm_filepath

        def load_from_store(gsm: GSM, self_: StoreElement):
            adata = read_h5ad(self_.get_path().joinpath("dataset.h5ad"))
            return adata

        self.store_files = [
            StoreElement(
                name=f"Wu_et_al_2021_raw_{gsm}",
                linux_hash=linux_hash_,
                darwin_hash=darwin_hash_,
                store=self.store,
                derivation=partial(derivation, gsm),
                derivation_files=[],
                load_from_store=partial(load_from_store, gsm)
            )
            for gsm, linux_hash_, darwin_hash_ in zip(self.gsm_list, self.linux_hashes, self.darwin_hashes)
        ]

    def get(self, gsm: GSM | None = None) -> AnnData:
        if isinstance(gsm, str):
            logging.debug("Turning string in GSM!")
            gsm = GSM(gsm)

        if gsm is None:
            adatas = []
            for store_file, gsm_ in zip(self.store_files, self.gsm_list):
                adata = store_file.get()
                assert "GSM" not in adata.obs.columns, "Dataset already has GSM!"
                adata.obs["GSM"] = repeat(str(gsm_), repeats=adata.n_obs)
                adatas.append(adata)
            adata = ad_concat(adatas, axis=0)
        elif gsm in self.gsm_list:
            logging.debug(f"GSM \"{gsm}\" is in gsm_list")
            adata = self.store_files[self.gsm_list.index(gsm)].get()
        else:
            raise ValueError("gsm must be one of GSM list or None!")
        return adata

    def exists(self) -> bool:
        return all([storeElement.exists() for storeElement in self.store_files])


def mutations_tcga(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        path = Path("../download/mutations.h5pd")
        copyfile(path, output_dir.joinpath("mutations.h5pd"))

    def load_from_store(self_: StoreElement) -> AnnData:
        return pd_read_hdf(self_.get_path().joinpath("mutations.h5pd"))

    store_element = StoreElement(
        name="mutations_TCGA",
        linux_hash=Hash("eb7b80a52bb3ba2840db66e3d5cef16e"),
        darwin_hash=Hash("eb7b80a52bb3ba2840db66e3d5cef16e"),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def tcga_samples(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        path = Path("../download/samples")
        for filename in listdir(path):
            adata = read_h5ad(path.joinpath(filename))
            number_of_genes = adata.n_obs
            adata = DataArray(
                adata.X.T,
                coords={"sample_barcode": adata.var.index, "Ensemble_Gene_ID": adata.obs.index}
            ).astype("int64")
            adata = adata.chunk({"sample_barcode": 1, "Ensemble_Gene_ID": number_of_genes})
            adata.to_netcdf(
                output_dir.joinpath(filename.removesuffix(".h5ad") + ".hdf5"),
                mode="w", format="NETCDF4", engine="netcdf4"
            )

    def load_from_store(self_: StoreElement) -> AnnData:
        data_arrays = {}
        for file in self_.get_path().iterdir():
            data_arrays[file.name.removesuffix(".hdf5")] = open_dataarray(file)
        return data_arrays

    store_element = StoreElement(
        name="TCGA_samples",
        linux_hash=Hash(None),
        darwin_hash=Hash("1974c5314603a1fe2703bb20879f9e6e"),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def deconvolution_fractions_tcga(store: Store):
    def derivation(output_dir, derivation_store_files):
        logger.debug("Deriving deconvolution fractions TCGA")
        assert len(derivation_store_files) == 0
        directory = Path("../download/references_ensembl")
        cell_state_to_cell_type_dict = {}
        for filepath in directory.iterdir():
            d = read_h5ad(filepath).uns["cell_state_to_cell_type_dict"]
            cell_state_to_cell_type_dict[filepath.name.removesuffix(".h5ad")] = {
                (v if isinstance(vs, ndarray) else str(vs)): k for k, vs in d.items() for v in vs}

        del cell_state_to_cell_type_dict["SKCM"]["{'Monocytes': 'Macrophages/Monocytes'}"]
        cell_state_to_cell_type_dict["SKCM"]["Macrophages/Monocytes"] = "Monocytes"

        directory = Path("../download/deconvolution_results/")
        for subdirectory in directory.iterdir():
            mkdir(output_dir.joinpath(subdirectory.name))
            for filepath in subdirectory.iterdir():
                if match(pattern=r"\.metadata\.md$", string=filepath.name):
                    continue
                adata = read_h5ad(filepath)
                data_array = DataArray(adata.X, dims=("cell_state", "gene_id"),
                                          coords={"cell_state": adata.obs_names, "gene_id": adata.var_names})
                gene_name_array = DataArray(adata.var["gene_name"], dims=("gene_id"),
                                               coords={"gene_id": adata.var_names})
                cell_type_array = DataArray(
                    [cell_state_to_cell_type_dict[subdirectory.name][cell_state] for cell_state in adata.obs_names],
                    dims=("cell_state"), coords={"cell_state": adata.obs_names})
                data = Dataset(
                    data_vars={
                        "expression_levels": data_array,
                        "gene_names": gene_name_array,
                        "cell_type": cell_type_array
                    }
                )
                data.to_netcdf(
                    output_dir.joinpath(subdirectory.name, filepath.name.removesuffix(".h5ad") + ".hdf5"),
                    mode="w", format="NETCDF4", engine="netcdf4"
                )
                logger.debug("Wrote deconvolution fractions TCGA for store")

    def load_from_store(self_: StoreElement) -> AnnData:
        logger.debug("Loading deconvolution fractions TCGA from store")
        data = {}
        for directory in self_.get_path().iterdir():
            logger.debug(f"Loading deconvolution fractions TCGA {directory.name} from store")
            data[directory.name] = {}
            for file in directory.iterdir():
                data[directory.name][file.name.removesuffix(".hdf5")] = open_dataset(file)
        logger.debug("Loaded deconvolution fractions TCGA from store")
        return data

    store_element = StoreElement(
        name="deconvolution_fractions_tcga",
        linux_hash=Hash("969f0a0040d2d1293e0f4f06b565f3fb"),
        darwin_hash=Hash("d1f36459248eeddb3797eb1772ef92de"),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def deconvolution_fractions_tcga_all(store: Store):
    def derivation(output_dir, derivation_store_files):
        logger.debug("Combining deconvolution fractions TCGA")
        assert len(derivation_store_files) == 1
        data = derivation_store_files[0].get()
        data_tree = {}
        for project_name, project_data in data.items():
            logger.debug(f"Combining deconvolution fractions {project_name}")
            expression_datas = []
            barcodes = []
            for k, v in project_data.items():
                barcodes.append(k)
                expression_datas.append(v)
            project_data = xr_concat(expression_datas, dim="barcode")
            project_data = project_data.assign_coords({"barcode": barcodes})
            project_data["gene_names"] = project_data["gene_names"].isel(barcode=0)
            project_data["cell_type"] = project_data["cell_type"].isel(barcode=0)
            data_tree[project_name] = DataTree(project_data)
        logger.debug(f"Combining deconvolution fractions into data tree")
        DataTree(children=data_tree).to_netcdf(
            output_dir.joinpath("deconvolution_fractions_tcga_all.hdf5"),
            mode="w", format="NETCDF4", engine="netcdf4"
        )

    def load_from_store(self_: StoreElement) -> AnnData:
        return open_datatree(self_.get_existing_path().joinpath("deconvolution_fractions_tcga_all.hdf5"), engine="netcdf4")

    store_element = StoreElement(
        name="deconvolution_fractions_tcga_all",
        linux_hash=Hash("d48e941d408ad73173a3f2539f6c68fb"),
        darwin_hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[deconvolution_fractions_tcga(store=store)],
        load_from_store=load_from_store
    )
    return store_element


def deconvolution_fractions_brca(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        path = Path("../download/deconvolution_results/BRCA")
        files = [file for file in listdir(path) if file.endswith(".h5ad")]
        assert len(files) == len(unique(files))

        fractions = []
        for file in files:
            adata = read_h5ad(path.joinpath(file))
            fraction_df = calculate_fractions(adata)
            fraction_df.rename(columns={"fraction": file[:-5]}, inplace=True)
            fractions.append(fraction_df.T)

        fractions = pd_concat(fractions)
        fractions.to_hdf(output_dir.joinpath("deconvolution_results.h5pd"), key="deconvolution_results")

    def load_from_store(self_: StoreElement) -> AnnData:
        return pd_read_hdf(self_.get_path().joinpath("deconvolution_results.h5pd"))

    store_element = StoreElement(
        name="deconvolution_fractions_brca",
        linux_hash=Hash(None),
        darwin_hash=Hash("233867510f5607e7c49e3cbbb29921ff"),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def mean_fractions_by_mutations_brca(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 2
        fractions = derivation_store_files[0].get()
        mutations = derivation_store_files[1].get()

        mutated_genes = unique(mutations["gene_ensemble_id"])

        data = Dataset(
            data_vars={
                "mean": DataArray(
                    dims=("gene", "cell_state"),
                    coords={
                        "gene": list(mutated_genes),
                        "cell_state": list(fractions.columns)
                    }),
                "std": DataArray(
                    dims=("gene", "cell_state"),
                    coords={
                        "gene": list(mutated_genes),
                        "cell_state": list(fractions.columns)
                    }),
                "number_of_samples": DataArray(
                    data=zeros((len(mutated_genes)), dtype=int),
                    dims="gene",
                    coords={
                        "gene": list(mutated_genes)
                    }
                )
            },
            coords={
                "gene": list(mutated_genes),
                "cell_state": list(fractions.columns)
            })

        mutation_to_barcodes = OrderedDict((mutated_gene, []) for mutated_gene in mutated_genes)
        for _, gene, barcode in mutations[["gene_ensemble_id", "barcode"]].itertuples():
            mutation_to_barcodes[gene].append(barcode)

        for mutated_gene in mutation_to_barcodes.keys():
            patients = patient(mutation_to_barcodes[mutated_gene])

            if len(patients) > 0:
                df = pd_concat([
                    fractions[array(patient(fractions.index)) == patient_] for patient_ in patients
                ])
                df_mean = df.mean(axis=0)
                df_std = df.std(axis=0)
                data["number_of_samples"].loc[mutated_gene] = df.shape[0]
                for cell_state in df.columns:
                    data["mean"].loc[mutated_gene, cell_state] = df_mean[cell_state]
                    data["std"].loc[mutated_gene, cell_state] = df_std[cell_state]

        # WARNING: engine="h5netcdf" gives nondeterministic files, ergo nondeterministic hashes!!!
        data.to_netcdf(output_dir.joinpath("mean_expressionsll.hdf5"), mode="w", format="NETCDF4", engine="netcdf4")

    def load_from_store(self_: StoreElement) -> AnnData:
        return open_dataset(self_.get_path().joinpath("mean_expressions.hdf5"))

    store_element = StoreElement(
        name="mean_fractions_by_mutations_brca",
        linux_hash=Hash(None),
        darwin_hash=Hash("4029303dbd2eee1b7ea96415b375dae2"),
        store=store,
        derivation=derivation,
        derivation_files=[deconvolution_fractions_brca(store), mutations_tcga(store)],
        load_from_store=load_from_store
    )
    return store_element


def wu_et_al_2021(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 1
        adata = derivation_store_files[0].get()
        adata.write_h5ad(output_dir.joinpath("dataset.h5ad"))

    def load_from_store(self_: StoreElement) -> AnnData:
        return read_h5ad(self_.get_path().joinpath("dataset.h5ad"))

    store_element = StoreElement(
        name="Wu_et_al_2021",
        linux_hash=Hash(None),
        darwin_hash=Hash("4f261fb51c4992aea7e14848fe3cac96"),
        store=store,
        derivation=derivation,
        derivation_files=[WuEtAl2021(store=store)],
        load_from_store=load_from_store
    )
    return store_element


def wu_et_al_2021_centroids(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 1
        adata = derivation_store_files[0].get()
        adata = adata[~adata.obs_names.duplicated()].copy()  # remove duplicates
        centroids = calc_centroids(adata, n_clusters=10, seed=0)
        centroids.write_h5ad(output_dir.joinpath("centroids.h5ad"))

    def load_from_store(self_: StoreElement) -> AnnData:
        return read_h5ad(self_.get_path().joinpath("centroids.h5ad"))

    store_element = StoreElement(
        name="WuEtAl2021_Centroids",
        hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[WuEtAl2021(store=store)],
        load_from_store=load_from_store
    )
    return store_element


def petralia_et_al_2024_raw(store: Store):
    """
    Patient data from CPTAC python package and supplementary data in downloads directory
    https://paynelab.github.io/cptac/tutorial01_data_intro.html
    :param store:
    :return:
    """

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0

    def load_from_store(self_: StoreElement) -> AnnData:
        return None

    store_element = StoreElement(
        name="Petralia_et_al_2024_raw",
        hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def thorsson_et_al_2018_raw(store: Store):
    """
    Data to cluster TMEs into immune subtypes
    :param store:
    :return:
    """

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0

    def load_from_store(self_: StoreElement) -> AnnData:
        return None

    store_element = StoreElement(
        name="Thorsson_et_al_2018_raw",
        hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def finotello_et_al_2019_raw(store: Store):
    """
    Signature for deconvolution
    :param store:
    :return:
    """

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0

    def load_from_store(self_: StoreElement) -> AnnData:
        return None

    store_element = StoreElement(
        name="Finotello_et_al_2019_raw",
        hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def bagaev_et_al_2021_signature_genes_for_tme_classification(store: Store):

    """
    Title: Conserved pan-cancer microenvironment subtypes predict response to immunotherapy

    MORE DATASETS!!! SEE FOLLOWING QUOTE:
    In two independent skin cutaneous melanoma cohorts (n = 58) (Nathanson et al., 2017; Snyder et al., 2014; Van Allen
    et al., 2015) treated with anti-CTLA-4 therapy, patients were classified into the four TME subtypes (Figure 5A), which
    significantly correlated with response to ipilimumab. The percentage of responders (R + long-term survivors [LS]) to
    anti-CTLA-4 therapy in the immune, non-fibrotic TME subtype IE was 82% in contrast to only 10% of subtype F (Figure 5B).
    In both cohorts, OS following antiCTLA-4 was the longest in TME subtype IE (Figures 5C and S8A). Similar findings were
    observed with three independent anti-CTLA-4 naive cohorts of melanoma patients (n = 114) treated with anti-PD1 therapy
    (Figure 5D) (Gide et al., 2019; Hugo et al., 2016; Liu et al., 2019a). Overall response was significantly higher in
    patients with TME subtype IE (75%) compared with subtype F (10%) (Figure 5E), with prolonged PFS and OS also noted in
    TME subtype IE in combined and individual cohort analysis (Figures 5F, S8B, and S8C)
    """

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        df = read_excel(
            "../download/bagaev_et_al_2021-conserved_pan-cancer_microenvironment/1-s2.0-S1535610821002221-mmc2.xlsx",
            sheet_name="Gene Description",
            skiprows=1
        )
        df.to_hdf(output_dir.joinpath("signature_gene_sets.h5pd"), key="signature_gene_sets")

    def load_from_store(self_: StoreElement) -> AnnData:
        return pd.read_hdf(self_.get_path().joinpath("signature_gene_sets.h5pd"))

    store_element = StoreElement(
        name="Bagaev_et_al_2021_signature_genes_for_tme_classification",
        linux_hash=Hash(None),
        darwin_hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def coulton_et_al_2024(store: Store):

    """
    Title: Using a pan-cancer atlas to investigate tumour associated macrophages as regulators of immunotherapy response
    Data from: https://zenodo.org/records/11222158
    """

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0

    def load_from_store(self_: StoreElement) -> AnnData:
        return None

    store_element = StoreElement(
        name="Coulton_et_al_2024",
        hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def eraslan_et_al_2022(store: Store):

    """
    Title: Single-nucleus cross-tissue molecular reference maps toward understanding disease gene function
    Data from: https://gtexportal.org/home/downloads/adult-gtex/single_cell
    Atlas of snRNAseq data of helthy tissue from all around the human body

    Example Code:
        filename_hdf = "../download/coulton_et_al_2024/tumor_associated_macrophages.h5Seurat"
        with h5py.File(filename_hdf, 'r') as hf:
            print(hf["active.ident"]["levels"], "| *Louvain cluster names*") # Louvain clusters
            print(hf["active.ident"]["values"], "| *Louvain cluster indexes of each cell???*") # Louvain clusters
            print(hf["assays"]["RNA"]["data"]["data"], "| *log RNA counts???*")
            print(hf["assays"]["RNA"]["data"]["indptr"], "| *some index???*")
            print(hf["assays"]["RNA"]["features"], "| *genes*")
            print(hf["assays"]["RNA"]["data"]["data"][:10])
            print(np.round((np.exp(hf["assays"]["RNA"]["data"]["data"][:10])-1)/4.08163265-1))
    """

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0

    def load_from_store(self_: StoreElement) -> AnnData:
        return None

    store_element = StoreElement(
        name="Eraslan_et_al_2022",
        hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def charoentong_et_al_2017(store: Store):
    """Signature genes for immune cell types (including subtypes)"""

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        df = read_excel("../download/Charoentong_et_al_2017/1-s2.0-S2211124716317090-mmc3.xlsx", skiprows=2)
        marker_genes_for_cell_types = OrderedDict()
        for _, gene, cell_type, _ in df.itertuples():
            if cell_type not in marker_genes_for_cell_types.keys():
                marker_genes_for_cell_types[cell_type] = []
            marker_genes_for_cell_types[cell_type].append(gene)

        with open(output_dir.joinpath("marker_genes_for_cell_types.pkl"), "wb") as f:
            dump(marker_genes_for_cell_types, f)

    def load_from_store(self_: StoreElement):
        with open(self_.get_path().joinpath("marker_genes_for_cell_types.pkl"), 'rb') as f:
            marker_genes_for_cell_types = load(f)
        return marker_genes_for_cell_types

    store_element = StoreElement(
        name="Charoentong_et_al_2017",
        linux_hash=Hash(None),
        darwin_hash=Hash("17bb9892379a670c7daac0acfc40e039"),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def charoentong_et_al_2024_interesting_genes(store: Store):
    """Signature genes for immune cell types (including subtypes)"""

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        interesting_genes = {
            "Antigen_processing": ["B2M", "TAP1", "TAP2", "HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DPB2",
                                   "HLA-DQA1", "HLA-DQB1", "HLA-DQB2", "HLA-DRB", "HLA-E", "HLA-F", "HLA-G"],
            "Immunoinhibitors": ["CTLA4", "PD1", "LAG3", "BTLA", "CD160", "IDO1", "IL10", "TIGIT", "TGFB1", "TIM3",
                                 "PDL1", "PDL2"],
            "Immunostimulators": ["CD27", "CD28", "C40LG", "CD40", "CD70", "CD80", "CD86", "ICOS", "IL6", "TMEM173",
                                  "TNFRSF13B", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFRSF13B",
                                  "TNFRSF13"]
        }

        with open(output_dir.joinpath("interesting_genes.pkl"), "wb") as f:
            dump(interesting_genes, f)

    def load_from_store(self_: StoreElement):
        with open(self_.get_path().joinpath("interesting_genes.pkl"), 'rb') as f:
            interesting_genes = load(f)
        return interesting_genes

    store_element = StoreElement(
        name="Charoentong_et_al_2017",
        linux_hash=Hash(None),
        darwin_hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def mlynska_et_al_2024_interesting_genes(store: Store):
    """
    Signature genes for immune cell types (including subtypes)

    TME classification matrix (see fig 2c):

            DESERT      |       EXCLUDED        |       INFLAMED        ||
            ====================================================================================
            low         |       mid             |       low             ||       ANGIOGENESIS
            very low    |       very low        |       very high       ||       IMMMUNE RESPONSE
            low         |       high            |       mid             ||       STROMA
    """


    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        interesting_genes = {
            "angiogenesis": ["ESM1", "FLT1", "KDR", "VEGFA"],
            "immune_response": ["CD8A", "CD27", "CD274", "CTLA4", "CXCL10", "CXCL11", "EOMES", "GZMA", "GZMB", "IDO1", "PDCD1",
                   "PRF1", "PSMB8", "PSMB9", "TAP1", "TAP2", "TIGIT"],
            "stroma": ["CXCL8", "CD163", "CD4", "CD68", "COL5A1", "COL5A2", "FAP", "FCGR2A", "FN1", "IL10", "IL1B", "IL6",
                   "LOXL1", "MMP9", "MS4A4A", "POSTN", "PTGS2", "TDO2"]
        }

        with open(output_dir.joinpath("interesting_genes.pkl"), "wb") as f:
            dump(interesting_genes, f)

    def load_from_store(self_: StoreElement):
        with open(self_.get_path().joinpath("interesting_genes.pkl"), 'rb') as f:
            interesting_genes = load(f)
        return interesting_genes

    store_element = StoreElement(
        name="mlynska_et_al_2024",
        linux_hash=Hash(None),
        darwin_hash=Hash("bcbd8d2d075485d92034fc8df49ad457"),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def ma_et_al_2022_macrophage_markers(store: Store):
    """Signature genes for different macrophage subtype"""

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0

        df = pd.read_html("/Users/halu/tmp/showFullTableHTML.html")[0]

        df["is_protein_marker"] = [
            False, True, False,  # INF
            False, False,  # Inflam
            False, True, False,  # LA
            False, True, False,  # Angio
            False, False, True,  # Reg
            False, False,  # Prolif
            False, True,  # RTM-TAMs - Kupffer cell-like
            False, False,  # RTM-TAMs - Alveolar RTM-like
            False, True, False,  # RTM-TAMs - Microglia-like
            False, True,  # RTM-TAMs - MT-RTM-like
            False,  # RTM-TAMs - Interstitial RTM-like
            False, True, False, True,  # Classical TIMs
            False, True, False, True,  # Nonclassical monocytes
            False, True  # Intermediate monocytes
        ]
        df["Protein_Marker"] = empty(df.shape[0], dtype=str)
        for protein_row in df[array(df["is_protein_marker"]) == True].iterrows():
            for normal_row in df.iterrows():
                if (protein_row[1].drop(["Signature", "is_protein_marker", "Protein_Marker"]) == normal_row[1].drop(
                        ["Signature", "is_protein_marker", "Protein_Marker"])).all():
                    df.loc[normal_row[0], "Protein_Marker"] = protein_row[1]["Signature"]

        df = df.drop(df[df["is_protein_marker"]].index)
        df = df.drop("is_protein_marker", axis=1)
        columns = df.columns
        df = df[concatenate((columns[:4], [columns[-1]], list(columns[4:-1])))]

        df.to_hdf(output_dir.joinpath("macrophage_markers.h5pd"))

    def load_from_store(self_: StoreElement):
        df = read_hdf(self_.get_path().joinpath("macrophage_markers.h5pd"))
        return df

    store_element = StoreElement(
        name="Charoentong_et_al_2017",
        linux_hash=Hash(None),
        darwin_hash=Hash("17bb9892379a670c7daac0acfc40e039"),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def sc_best_practices_cell_markers(store: Store):
    """https://www.sc-best-practices.org/cellular_structure/annotation.html"""

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        marker_genes_coarse = {
            "CD14+ Mono": ["FCN1", "CD14"],
            "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
            # Note: DMXL2 should be negative
            "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
            "Erythroblast": ["MKI67", "HBA1", "HBB"],
            # Note HBM and GYPA are negative markers
            "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
            "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
            "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
            "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
            # Note IGHD and IGHM are negative markers
            "B cells": [
                "MS4A1",
                "ITGB1",
                "COL4A4",
                "PRDM1",
                "IRF4",
                "PAX5",
                "BCL11A",
                "BLK",
                "IGHD",
                "IGHM",
            ],
            "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
            # Note PAX5 is a negative marker
            "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
            "CD4+ T": ["CD4", "IL7R", "TRBC2"],
            "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
            "T naive": ["LEF1", "CCR7", "TCF7"],
            "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
        }
        marker_genes_fine = {
            "CD14+ Mono": ["FCN1", "CD14"],
            "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
            "ID2-hi myeloid prog": [
                "CD14",
                "ID2",
                "VCAN",
                "S100A9",
                "CLEC12A",
                "KLF4",
                "PLAUR",
            ],
            "cDC1": ["CLEC9A", "CADM1"],
            "cDC2": [
                "CST3",
                "COTL1",
                "LYZ",
                "DMXL2",
                "CLEC10A",
                "FCER1A",
            ],  # Note: DMXL2 should be negative
            "Normoblast": ["SLC4A1", "SLC25A37", "HBB", "HBA2", "HBA1", "TFRC"],
            "Erythroblast": ["MKI67", "HBA1", "HBB"],
            "Proerythroblast": [
                "CDK6",
                "SYNGR1",
                "HBM",
                "GYPA",
            ],  # Note HBM and GYPA are negative markers
            "NK": ["GNLY", "NKG7", "CD247", "GRIK4", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
            "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
            "Lymph prog": [
                "VPREB1",
                "MME",
                "EBF1",
                "SSBP2",
                "BACH2",
                "CD79B",
                "IGHM",
                "PAX5",
                "PRKCE",
                "DNTT",
                "IGLL1",
            ],
            "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
            "B1 B": [
                "MS4A1",
                "SSPN",
                "ITGB1",
                "EPHA4",
                "COL4A4",
                "PRDM1",
                "IRF4",
                "CD38",
                "XBP1",
                "PAX5",
                "BCL11A",
                "BLK",
                "IGHD",
                "IGHM",
                "ZNF215",
            ],  # Note IGHD and IGHM are negative markers
            "Transitional B": ["MME", "CD38", "CD24", "ACSM3", "MSI2"],
            "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
            "Plasmablast": ["XBP1", "RF4", "PRDM1", "PAX5"],  # Note PAX5 is a negative marker
            "CD4+ T activated": ["CD4", "IL7R", "TRBC2", "ITGB1"],
            "CD4+ T naive": ["CD4", "IL7R", "TRBC2", "CCR7"],
            "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
            "T activation": ["CD69", "CD38"],  # CD69 much better marker!
            "T naive": ["LEF1", "CCR7", "TCF7"],
            "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
            "G/M prog": ["MPO", "BCL2", "KCNQ5", "CSF3R"],
            "HSC": ["NRIP1", "MECOM", "PROM1", "NKAIN2", "CD34"],
            "MK/E prog": [
                "ZNF385D",
                "ITGA2B",
                "RYR3",
                "PLCB1",
            ],  # Note PLCB1 is a negative marker
        }
        marker = {"coarse": marker_genes_coarse, "fine": marker_genes_fine}
        with open(output_dir.joinpath("cell_markers.pickle"), 'wb') as f:
            dump(marker, f)

    def load_from_store(self_: StoreElement):
        with open(self_.get_path().joinpath("cell_markers.pickle"), 'rb') as f:
            markers = load(f)
        return markers

    store_element = StoreElement(
        name="SC-Best-Practices_Cell_Markers",
        linux_hash=Hash(None),
        darwin_hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element


def marija_hotspots(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 0
        directory = Path("../download/marija-hotspots/")

        genes_of_interest = pd.read_csv(directory.joinpath("genes_of_interest.txt"), header=None)
        genes_of_interest = genes_of_interest.rename(columns={0: "gene"})
        genes_of_interest.to_hdf(output_dir.joinpath("genes_of_interest.h5pd"), key="genes_of_interest")

        genes_of_special_interest = pd.read_csv(directory.joinpath("genes_of_special_interest.txt"), header=None)
        genes_of_special_interest = genes_of_special_interest.rename(columns={0: "gene"})
        genes_of_special_interest.to_hdf(output_dir.joinpath("genes_of_special_interest.h5pd"), key="genes_of_special_interest")

        hotspot_genes = pd.read_csv(directory.joinpath("hotspot_genes.tsv"), sep="\t")
        hotspot_genes.to_hdf(output_dir.joinpath("hotspot_genes.h5pd"), key="hotspot_genes")

        publication_hotspots = pd.read_csv(directory.joinpath("publication_hotspots.tsv"), sep="\t")
        publication_hotspots.to_hdf(output_dir.joinpath("publication_hotspots.h5pd"), key="publication_hotspots")


    def load_from_store(self_: StoreElement):
        return {
            key: pd.read_hdf(self_.get_path().joinpath(key + ".h5pd"))
            for key in ("genes_of_interest", "genes_of_special_interest", "hotspot_genes", "publication_hotspots")
        }

    store_element = StoreElement(
        name="Marija_Hotspot-Genes",
        linux_hash=Hash(None),
        darwin_hash=Hash(None),
        store=store,
        derivation=derivation,
        derivation_files=[],
        load_from_store=load_from_store
    )
    return store_element
