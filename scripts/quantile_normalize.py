import pandas as pd
import numpy as np
from utility.ryp import r, to_r, to_py
from src.quantile_normalize import quantile_normalize

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

def get_results(RE, algorithm):
    """Get results from raw RE scores

    Args:
        RE: pd.DataFrame
        algorithm: str, either 'R', 'nxu' or 'huosan'

    Returns:
        pd.DataFrame
    """
    if algorithm == 'R':
        RE_normalized = qn_R(RE)
    elif algorithm == 'nxu':
        RE_normalized = qn_python(RE)
    elif algorithm == 'huosan':
        RE_normalized = pd.DataFrame(quantile_normalize(RE),index = RE.index, columns = RE.columns)

    # get results
        
    to_r(RE_normalized, "RE")
    r("""
    RE_KO <- RE[, stringr::str_detect(colnames(RE), "4_7|4_9|95_3")]
    RE_WT <- RE[, stringr::str_detect(colnames(RE), "87|90|91")]

    log2_RE_FC <- log(rowMeans(RE_KO, na.rm = TRUE) / rowMeans(RE_WT, na.rm = TRUE), 2)
    Z_score <- rowMeans(RE_KO, na.rm = TRUE) - rowMeans(RE_WT, na.rm = TRUE) / sqrt(apply(RE_KO, 1, sd, na.rm = TRUE)**2 + apply(RE_WT, 1, sd, na.rm = TRUE)**2)

    results <- data.frame(gene_id = row.names(RE), log2_RE_FC = log2_RE_FC, Z_score = Z_score)
      """
    )    
    results = to_py("results", format='pandas')        
    return results

# load data

RE = pd.read_csv("results/RE_quant/linux_RE_unormalized.csv", index_col=0)
RE_clean = remove_inf_nan(RE)

pd.DataFrame(quantile_normalize(RE),index = RE.index, columns = RE.columns)

get_results(RE, 'R')
get_results(RE, 'huosan')

to_py("RE_KO", format='pandas')
to_py("RE_WT", format='pandas')
to_py("log2_RE_FC", format='pandas')
r('rowMeans(RE_KO, na.rm = TRUE) %>% head()')
r('rowMeans(RE_WT, na.rm = TRUE) %>% head()')
r('log(rowMeans(RE_KO, na.rm = TRUE) / rowMeans(RE_WT, na.rm = TRUE), 2) %>% head()')

# Withiout cleaning

RE_normalized_py = qn_python(RE)
RE_normalized_R = qn_R(RE)

# With cleaning

RE_normalized_py = qn_python(RE_clean)
RE_normalized_R = qn_R(RE_clean)

# Plotting
import seaborn as sns
import matplotlib.pyplot as plt

pd.melt(RE_normalized_py)\
    .assign(value = lambda x: np.log2(x['value']+1))\
    .pipe((sns.boxplot, 'data'), x='variable', y='value')\
    .set_title('Quantile normalized expression values')