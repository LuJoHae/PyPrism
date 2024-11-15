import logging
from pathlib import Path
import requests
import GEOparse
from io import TextIOWrapper
from tempfile import TemporaryDirectory
import os
import hashlib
import pickle
import sys
from typing import List
from checksumdir import dirhash
import re
import shutil
import tarfile
from anndata import read_mtx
from pandas import read_csv
import weakref
from sys import platform


class WeakMethod:
    def __init__(self, func, instance):
        self.func = func
        self.instance_ref = weakref.ref(instance)

        self.__wrapped__ = func  # this makes things like `inspect.signature` work

    def __call__(self, *args, **kwargs):
        instance = self.instance_ref()
        return self.func(instance, *args, **kwargs)

    def __repr__(self):
        cls_name = type(self).__name__
        return '{}({!r}, {!r})'.format(cls_name, self.func, self.instance_ref())


class HiddenConsoleOutput:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        self._original_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        sys.stderr.close()
        sys.stderr = self._original_stderr


class Store:
    def __init__(self, path: Path = Path("store")):
        self._path = Path(path)

    def __str__(self):
        return f"Store at \"{self._path}\""

    def get_path(self):
        return self._path

    def exists(self):
        return self._path.is_dir()

    def create(self):
        if self.exists():
            raise FileExistsError(f"Store already exists at \"{self.get_path()}\"")
        os.makedirs(self.get_path(), exist_ok=False)


class Hash:
    def __init__(self, hash: str):
        self._hash = hash

    def __str__(self):
        return str(self._hash)

    def __eq__(self, other):
        return self._hash == other._hash


class OSDependentHash:
    def __init__(self, linux_hash: Hash, darwin_hash: Hash):
        assert isinstance(linux_hash, Hash)
        assert isinstance(darwin_hash, Hash)
        self._linux_hash = linux_hash
        self._darwin_hash = darwin_hash

    def get_hash(self) -> Hash:
        if platform == "linux" or platform == "linux2":
            return self._linux_hash
        elif platform == "darwin":
            return self._darwin_hash
        else:
            raise OSError(f"Unsupported platform: {platform}")

    def set_hash(self, hash: Hash):
        if platform == "linux" or platform == "linux2":
            self._linux_hash = hash
        elif platform == "darwin":
            self._darwin_hash = hash
        else:
            raise OSError(f"Unsupported platform: {platform}")


class StoreElement:
    def __init__(self, name: str, linux_hash: Hash, darwin_hash: Hash, store: Store, derivation, derivation_files, load_from_store):
        self._name = name
        self._os_dependent_hash = OSDependentHash(linux_hash=linux_hash, darwin_hash=darwin_hash)
        self._store = store
        self._derivation = derivation
        self._derivation_store_files = derivation_files
        self._load_from_store = WeakMethod(load_from_store, self)

    def __str__(self):
        return f"File \"{self._name}\" in store \"{self._store}\" with hash \"{str(self.get_hash())}\""

    def get_hash(self) -> Hash:
        return self._os_dependent_hash.get_hash()

    def set_hash(self):
        return self._os_dependent_hash.set_hash()

    def get_path(self):
        return self._store.get_path() / f"{str(self.get_hash())}-{self._name}"

    def store_exists(self):
        return self._store.exists()

    def exists(self):
        if not self.store_exists():
            raise FileNotFoundError("Store does not exist!")
        return self.get_path().exists()

    def _download(self):
        with TemporaryDirectory() as tmpdir:
            logging.debug(f"Temporary directory: {tmpdir}")
            filepath = self._derivation(output_dir=Path(tmpdir), derivation_store_files=self._derivation_store_files)
            true_hash = Hash(dirhash(tmpdir))
            if self.get_hash()._hash is None:
                logging.debug("Hash is None and hash checking is ignored!")
                self.set_hash(true_hash)
            assert self.get_hash() == true_hash, f"Hash is incorrect! Expected hash: {self.get_hash()}, got {true_hash}"
            logging.debug(f"Hash: {self.get_hash()}")
            os.replace(tmpdir, self.get_path())

    def _download_derivation_files(self):
        logging.debug(f"{self} does not exist and must be derived!")
        for derivation_store_file in self._derivation_store_files:
            logging.debug(f"Derivation {derivation_store_file} does not exist and must be derived!")
            if not derivation_store_file.exists():
                derivation_store_file.get()

    def get(self):
        if self.exists():
            logging.debug(f"Store file {self} already exists!")
        else:
            self._download_derivation_files()
            self._download()
        logging.debug(f"load_from_store: {self._load_from_store}")
        return self._load_from_store()


class StoreElementContainer:
    def __init__(self):
        pass


class URL:
    def __init__(self, url: str):
        self._url = url

    def __str__(self):
        return str(self._url)


class GSM:
    def __init__(self, gsm: str):
        self._gsm = gsm

    def __str__(self):
        return self._gsm

    def __eq__(self, other):
        return self._gsm == other._gsm


def get_gsm_tar_filepath(directory):
    subdirs = os.listdir(directory)
    assert len(subdirs) == 1, "More that one directory or file in directory!"
    subdir = directory.joinpath(subdirs[0])
    assert os.path.isdir(subdir), "Subdir is not a directory!"
    files = os.listdir(subdir)
    assert len(files) == 1, "More that one file in subdir!"
    filepath = subdir.joinpath(files[0])
    assert os.path.isfile(filepath), "File is not a file!"
    assert filepath.suffixes == [".tar", ".gz"], "File extension is not .tar.gz!"
    return filepath


def get_path_singleton(path: Path):
    files = os.listdir(path)
    assert len(files) == 1, "More that one file in path!"
    return path.joinpath(files[0])


def extract_and_delete_directory(dir: Path):
    for file_name in os.listdir(dir):
        shutil.move(os.path.join(dir, file_name), dir.parent)
    Path.rmdir(dir)


def download_gsm(gsm: GSM, output_dir: TextIOWrapper) :
    output_dir = Path(output_dir)
    logging.debug(f"output_dir: {output_dir}")
    gsm = GEOparse.get_GEO(geo=gsm._gsm, destdir="output_dir", silent=True)
    with HiddenConsoleOutput():
        gsm.download_supplementary_files(directory=output_dir, download_sra=False)
    supplementary_directory = get_path_singleton(output_dir)
    logging.debug(f"supplementary_directory.name: {supplementary_directory.name}")
    assert re.compile("^Supp_GSM[0-9]{7}_CID[0-9]{4}[1,A,N]?$").match(supplementary_directory.name), \
        "Supplementary directory does not match regex!"
    extract_and_delete_directory(supplementary_directory)
    supplementary_tarfile = get_path_singleton(output_dir)
    logging.debug(f"supplementary_directory.name: {supplementary_directory.name}; supplementary_tarfile.name: {supplementary_tarfile.name}")
    assert supplementary_directory.name[5:] + ".tar.gz" == supplementary_tarfile.name
    with tarfile.open(supplementary_tarfile, "r") as tar:
        tar.extractall(output_dir)
    os.remove(supplementary_tarfile)
    cid_directory = get_path_singleton(output_dir)
    logging.debug(f"supplementary_directory.name[17:]: {supplementary_directory.name[16:]}; cid_directory.name: {cid_directory.name}")
    assert cid_directory.name == supplementary_directory.name[16:]
    extract_and_delete_directory(cid_directory)

    count_filepath = output_dir.joinpath("count_matrix_sparse.mtx")
    adata_filepath = output_dir.joinpath("dataset.h5ad")
    assert Path.exists(output_dir.joinpath(count_filepath))
    adata = read_mtx(count_filepath).T
    barcodes = read_csv(output_dir.joinpath("count_matrix_barcodes.tsv"), sep="\t", header=None)
    genes = read_csv(output_dir.joinpath("count_matrix_genes.tsv"), sep="\t", header=None)
    logging.debug(f"genes.columns: {genes.columns}")
    logging.debug(genes)
    metadata = read_csv(output_dir.joinpath("metadata.csv"), sep=",", header=0)
    assert all(metadata.iloc[:, 0] == barcodes.iloc[:, 0]), "Barcodes in metadata are incorrect!"
    metadata.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
    metadata.set_index("barcode", inplace=True)
    logging.debug(f"metadata.columns: {metadata.columns}")
    adata.obs = metadata
    adata.var_names = genes[0]
    adata.uns["gsm_metadata"] = gsm.metadata
    adata.write_h5ad(adata_filepath, compression="gzip", compression_opts=4)
    os.remove(output_dir.joinpath("count_matrix_sparse.mtx"))
    os.remove(output_dir.joinpath("count_matrix_barcodes.tsv"))
    os.remove(output_dir.joinpath("count_matrix_genes.tsv"))
    os.remove(output_dir.joinpath("metadata.csv"))
    logging.debug(f"Directory list (at end): {os.listdir(output_dir)}")
