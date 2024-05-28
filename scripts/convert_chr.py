import pandas as pd
from pathlib import Path
import sys

junc_flist = [junc for junc in Path("proc/leafcutter/junc_files").iterdir()]

for junc_file in junc_flist:
    data = pd.read_csv(
        junc_file, 
        sep = "\t", 
        names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand"], 
        dtype = {"chrom": str},
        low_memory=False)

    data = data[~data["chrom"].str.contains("\.")]
    data["chrom"] = "chr" + data["chrom"]
    data.to_csv(junc_file, sep = "\t", index = False, header = False)