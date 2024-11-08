import requests
import tarfile
import os
import hashlib
import logging
import anndata as ad
import pandas as pd
from dataclasses import dataclass


def file_hash_md5(filename):
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def download_file_raw(url, filepath):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(filepath, 'wb') as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)


def download_file(url, filename, check_hash_md5, directory="."):
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.isdir(directory):
        raise NotADirectoryError("Output directory is not a directory!")

    filepath = os.path.join(directory, filename)
    if os.path.exists(filepath):
        if file_hash_md5(filepath) != check_hash_md5:
            raise FileExistsError("File already exists but md5 hash does not fit!")
        logging.info("File already there and md5 hash matched!")
        return

    download_file_raw(url=url, filepath=filepath)
    return


def extract_tarfile(filepath, extract_directory):
    if not tarfile.is_tarfile(filepath):
        raise tarfile.ReadError("File is not a tarball!")

    with tarfile.open(filepath, "r") as tf:
        tf.extractall(path=extract_directory)


@dataclass
class DataSourceDetails:
    name: str
    url: str
    html_url: str
    hash_md5: str
    adata_hash: str
    convert_to_h5ad: callable
    read_h5ad: callable


def _breast_cancer_single_cell_convert_to_h5ad(directory):
    adata = ad.read_mtx(directory + "/data/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx").T
    barcodes = pd.read_csv(directory + "/data/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv", index_col=False,
                           header=None)
    genes = pd.read_csv(directory + "/data/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv", index_col=False,
                        header=None)
    adata.var_names = genes[0]
    #adata.obs = barcodes.apply(lambda x: pd.Series(x[0].split("_")), axis=1).rename(columns={0: "sample", 1: "barcode"})
    #adata.obs = pd.concat([adata.obs, metadata], axis=1)
    #adata.obs_names = barcodes[0]
    adata.obs = pd.read_csv(directory + "/data/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
    adata.obs.set_index("Unnamed: 0", inplace=True)
    adata.obs.index.rename("barcode", inplace=True)
    adata.write_h5ad(directory + "/data.h5ad", compression="gzip")


def get_dataset_source_details(dataset_name):
    dataset_source_details = {
        "breast_cancer_bulk": DataSourceDetails(
            name="breast_cancer_bulk",
            url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FbulkRNAseq%5Fraw%5Fcounts%2Etxt%2Egz",
            html_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
            hash_md5="887d02e472b44eec76d42c5d7fcd9f7a",
            adata_hash="bad0eb54201527d0dc7e41bbafd9990c",
            convert_to_h5ad=lambda directory: ad.read_csv(
                "download/887d02e472b44eec76d42c5d7fcd9f7a/data/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt",
                delimiter='\t'
            ).T.write_h5ad(filename="download/887d02e472b44eec76d42c5d7fcd9f7a/data.h5ad", compression="gzip"),
            read_h5ad=lambda: ad.read_h5ad("download/887d02e472b44eec76d42c5d7fcd9f7a/data.h5ad")
        ),
        "breast_cancer_single_cell": DataSourceDetails(
            name="breast_cancer_single_cell",
            url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file&file=GSE176078%5FWu%5Fetal%5F2021%5FBRCA%5FscRNASeq%2Etar%2Egz",
            html_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078",
            hash_md5="e6d9df4a627bdd60cb4898be2ede77cd",
            adata_hash="65a44fc230c25ae58dbad2b478a3327f",
            convert_to_h5ad=_breast_cancer_single_cell_convert_to_h5ad,
            read_h5ad=lambda: ad.read_h5ad("download/e6d9df4a627bdd60cb4898be2ede77cd/data.h5ad")
        )
    }
    return dataset_source_details[dataset_name]


def is_adata(directory, hash):
    filepath = os.path.join(directory, "data.h5ad")
    if os.path.exists(filepath):
        return file_hash_md5(os.path.join(directory, "data.h5ad")) == hash
    else:
        return False


def get_dataset(dataset_name):
    data_source_details = get_dataset_source_details(dataset_name=dataset_name)
    logging.info("Getting dataset: {}".format(data_source_details.name))
    filename = "data.tar"
    directory = os.path.join("download", data_source_details.hash_md5)

    if os.path.exists(os.path.join(directory, "data.h5ad")):
        if is_adata(directory=directory, hash=data_source_details.adata_hash):
            logging.info("adata already there and hash matched")
            return data_source_details.read_h5ad()
        else:
            raise FileExistsError("data.h5ad already ther, but hash does not match!")

    filepath = os.path.join(directory, filename)
    logging.info("Download file")
    download_file(url=data_source_details.url, filename=filename, check_hash_md5=data_source_details.hash_md5,
                  directory=directory)
    logging.info("Extract")
    extraction_directory = os.path.join(directory, "data")
    extract_tarfile(filepath=filepath, extract_directory=extraction_directory)
    logging.info("Convert to h5ad")
    data_source_details.convert_to_h5ad(directory)
    logging.info("Read h5ad")
    return data_source_details.read_h5ad()