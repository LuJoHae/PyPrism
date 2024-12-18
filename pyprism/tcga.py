from dataclasses import dataclass
from typing import List
from enum import Enum
from collections.abc import Iterable


TCGAProjectName = Enum("TCGAProjectName", [
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
    "UCEC", "UCS", "UVM"
])


DiseaseType = Enum("DiseaseType", [
    "Adenomas and Adenocarcinomas",
    "Adnexal and Skin Appendage Neoplasms",
    "Basal Cell Neoplasms",
    "Complex Epithelial Neoplasms",
    "Cystic, Mucinous and Serous Neoplasms",
    "Ductal and Lobular Neoplasms",
    "Epithelial Neoplasms, NOS",
    "Fibroepithelial Neoplasms",
    "Squamous Cell Neoplasms"
])


PrimarySite = Enum("PrimarySite", [
    "Bones, joints and articular cartilage of limbs",
    "Colon",
    "Connective, subcutaneous and other soft tissues",
    "Corpus uteri",
    "Kidney",
    "Meninges",
    "Other and unspecified male genital organs",
    "Other and unspecified parts of tongue",
    "Ovary",
    "Peripheral nerves and autonomic nervous system",
    "Retroperitoneum and peritoneum",
    "Stomach",
    "Uterus, NOS"
])


@dataclass(frozen=True)
class TCGAProject:
    name: str
    disease_type: List[str]
    primary_site: List[str]
    number_of_cases: int

acc = TCGAProject(name="ACC", disease_type=["Adenomas and Adenocarcinomas"], primary_site=["Adrenal gland"], number_of_cases=92)

@dataclass(frozen=True)
class TCGAProjectSingleCellSource():
    tcga_project: TCGAProject
    studies: List[str]


TCGA = TCGAProjectSingleCellSource

tcga_project_single_cell_sources = [
    TCGA("BRCA", []),
    TCGA("GBM", []),
    TCGA("OV", []),
    TCGA("LUAD", []),
    TCGA("UCEC", []),
    TCGA("KIRC", []),
    TCGA("HNSC", []),
    TCGA("LGG", []),
    TCGA("THCA", []),
    TCGA("LUSC", []),
    TCGA("PRAD", []),
    TCGA("SKCM", []),
    TCGA("COAD", []),
    TCGA("STAD", []),
    TCGA("BLCA", []),
    TCGA("LIHC", []),
    TCGA("CESC", []),
    TCGA("KIRP", []),
    TCGA("TGCT", []),
    TCGA("SARC", []),
    TCGA("LAML", []),
    TCGA("PAAD", []),
    TCGA("ESCA", []),
    TCGA("PCPG", []),
    TCGA("READ", []),
    TCGA("THYM", []),
    TCGA("KICH", []),
    TCGA("ACC", []),
    TCGA("MESO", []),
    TCGA("UVM", []),
    TCGA("DLBC", []),
    TCGA("UCS", []),
    TCGA("CHOL", [])
]


def patient(barcode: str | Iterable[str]) -> str | List[str]:
    def take12(barcode_: str) -> str:
        return barcode_[:12]

    if isinstance(barcode, str):
        return take12(barcode[:12])
    if isinstance(barcode, Iterable):
        return list(map(take12, barcode))
