import logging
from pyprism.store import (Store, StoreElement, StoreElementContainer, Hash, download_gsm, GSM)
from io import TextIOWrapper
from pathlib import Path
from functools import partial
from anndata import AnnData, concat, read_h5ad
from numpy import repeat
import scanpy as sc


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
            adata = concat(adatas, axis=0)
        elif gsm in self.gsm_list:
            logging.debug(f"GSM \"{gsm}\" is in gsm_list")
            adata = self.store_files[self.gsm_list.index(gsm)].get()
        else:
            raise ValueError("gsm must be one of GSM list or None!")
        return adata

    def exists(self) -> bool:
        return all([storeElement.exists() for storeElement in self.store_files])


def wu_et_al_2021(store: Store):

    def derivation(output_dir, derivation_store_files):
        assert len(derivation_store_files) == 1
        adata = derivation_store_files[0].get()
        adata.write_h5ad(output_dir.joinpath("dataset.h5ad"))

    def load_from_store(self_: StoreElement) -> AnnData:
        return read_h5ad(self_.get_path().joinpath("dataset.h5ad"))

    store_element = StoreElement(
        name="Wu_et_al_2021",
        hash=Hash("4f261fb51c4992aea7e14848fe3cac96"),
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