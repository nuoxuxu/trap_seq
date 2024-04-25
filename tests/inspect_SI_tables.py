# read excel as pandas dataframe
import pandas as pd
import numpy as np

salmon = pd.read_csv("salmon/salmon_results/4_7_Input_S10_quant/quant.sf", sep="\t")
salmon_2 = pd.read_csv("salmon_2/4_7_Input_S10_L004_R_quant.sf", sep="\t")

salmon = salmon["Name"].str.split(".", expand = True).iloc[:, 0]
salmon_2 = salmon_2["Name"].str.split("|", expand=True).iloc[:, 0].str.split(".", expand=True).iloc[:, 0]

salmon_2.isin(salmon).sum()

np.intersect1d(salmon.values, salmon_2.values).shape