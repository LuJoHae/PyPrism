import pandas as pd
import math


if __name__ == "__main__":
    df = pd.read_excel("../../mmc2.xlsx", sheet_name="SolidBlood Cancer (GMT)", skiprows=1)

    last_id = None
    for index in range(df.shape[0]):
        if str(df["ID"][index]) == "nan":
            #print(df.iat[index, 0])
            df.at[index, "ID"] = last_id
        else:
            last_id = df["ID"][index]

    gene_lists = []
    for genes_row in df.iloc[:, 2:].iterrows():
        gene_lists = [gene for gene in genes_row[1] if str(gene) != 'nan']
