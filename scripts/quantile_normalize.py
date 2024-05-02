import pandas as pd
import numpy as np
from utility.ryp import r, to_r, to_py
r("library(dplyr)")

def qn_R(RE):
    to_r(RE, "RE", format='matrix')
    r('RE_normalized <- preprocessCore::normalize.quantiles(RE, keep.names=TRUE)')
    RE = to_py("RE_normalized", format='pandas')
    return RE

def remove_inf_nan(RE):
    # remove rows that contain any missing values
    RE = RE.loc[RE.isna().sum(axis=1)==0]

    # remove rows that contain any infinte values
    RE = RE[RE.isin([np.inf, -np.inf]).sum(axis=1)==0]

    # remove genes that are not expressed in any sample
    # RE = RE.loc[RE.sum(axis=1)>0]
    return RE

def qn_python(RE):
    RE_sorted = pd.DataFrame()
    for column in RE.columns:
        RE_sorted[column] = RE[column].sort_values(ignore_index=True).values
    ranked_means = RE_sorted.mean(axis=1)

    RE_normalized = pd.DataFrame()
    for column in RE.columns:
        RE_normalized[column] = ranked_means.iloc[RE[column].rank(axis=0).values - 1].values
    RE_normalized.index = RE.index 
    return RE_normalized

# Test in Linux

## With cleaning

RE = pd.read_csv("results/RE_quant/linux_RE_unormalized.csv", index_col=0)
RE_clean = remove_inf_nan(RE)
RE_normalized_py = qn_python(RE_clean)
RE_normalized_R = qn_R(RE_clean)

## Withiout cleaning

RE_normalized_py = qn_python(RE)
RE_normalized_R = qn_R(RE)

to_r(RE_normalized_R, "RE")

r(
    """
RE_KO <- RE[, stringr::str_detect(colnames(RE), "4_7|4_9|95_3")]
RE_WT <- RE[, stringr::str_detect(colnames(RE), "87|90|91")]

log2_RE_FC <- log(rowMeans(RE_KO, na.rm = TRUE) / rowMeans(RE_WT, na.rm = TRUE), 2)
Z_score <- rowMeans(RE_KO, na.rm = TRUE) - rowMeans(RE_WT, na.rm = TRUE) / sqrt(apply(RE_KO, 1, sd, na.rm = TRUE)**2 + apply(RE_WT, 1, sd, na.rm = TRUE)**2)

results <- data.frame(gene_id = row.names(RE), log2_RE_FC = log2_RE_FC, Z_score = Z_score)
    """
)

r(
    """    
results %>%
    filter_if(is.numeric, all_vars(!is.na(.) & !is.infinite(.))) %>%
    filter(abs(Z_score) > 2 & abs(log2_RE_FC) > 0.2) %>%
    arrange(desc(Z_score)) %>%
    head()
    """
)

results = to_py("results", format='pandas')

remove_inf_nan(results)\
    .query('abs(Z_score) > 2 & abs(log2_RE_FC) > 0.2')

(RE_normalized.sum(axis=1) == 0).sum()


# Test in Windows
RE_clean.to_csv("proc/RE_clean.csv")

# Plotting
import seaborn as sns
import matplotlib.pyplot as plt

pd.melt(RE_normalized_py)\
    .assign(value = lambda x: np.log2(x['value']+1))\
    .pipe((sns.boxplot, 'data'), x='variable', y='value')\
    .set_title('Quantile normalized expression values')