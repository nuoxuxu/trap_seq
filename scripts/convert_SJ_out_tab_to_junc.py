import pandas as pd
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description = "Convert STAR SJ.out.tab to leafcutter junc files")
parser.add_argument("--STAR_res", type = str, help = "Path to STAR results")
parser.add_argument("--juncfiles", type = str, help = "Path to juncfiles.txt")
parser.add_argument("--groups_file", type = str, help = "Path to groups_file.txt")
args = parser.parse_args()

def main():
    for tab_file in [path for path in Path(args.STAR_res).glob("**/SJ.out.tab")]:
        tab = pd.read_csv(
            tab_file,
            sep = "\t",
            names = ["chromosome", "start", "end", "strand", "motif", "annotated", "unique_reads", "multi_reads", "max_overhang"],
            dtype = {"chromosome": str})

        tab = tab\
            .assign(name = range(1, len(tab) + 1))\
            .rename(columns = {"chromosome": "chrom", "start": "chromStart", "end": "chromEnd", "unique_reads": "score"})\
            .loc[:, ["chrom", "chromStart", "chromEnd", "name", "score", "strand"]]

        tab["strand"] = tab["strand"].replace({1: "+", 2: "-"})
        prefix = tab_file.parent.name
        tab.to_csv(f"proc/leafcutter/junc_files/{prefix}.junc", sep = "\t", index = False, header = False)

    with open(args.juncfiles, "w") as f:
        f.write("\n".join([str(path.absolute()) for path in Path("proc/leafcutter/junc_files").iterdir()]))

    bam_list = [path.parent.name for path in Path(args.STAR_res).glob("**/Aligned.sortedByCoord.out.bam")]

    pd.DataFrame(bam_list,  columns = ["bam_file"])\
        .assign(condition = lambda x: x["bam_file"].str.count("_"))\
        .assign(condition = lambda x: x["condition"].replace({3: "KO", 2: "WT"}))\
        .to_csv(args.groups_file, sep = " ", index = False, header = False)

if __name__ == "__main__":
    main()