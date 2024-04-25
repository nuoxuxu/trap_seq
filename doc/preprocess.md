# Preprocessing
On SciNet Niagara clusters
1. Create `Homo_sapiens.GRCh38.cdna.ncrna.fa` by merging two fasta files. From [salmon_kalliso_STAR_compare.md](https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/salmon_kalliso_STAR_compare.md).
```bash
wget ftp://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
## merge together 
gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.cdna.ncrna.fa
```

2. Add eGFP nucleotide sequence to the FASTA file obtained from [EGFP SnapGene](https://www.snapgene.com/plasmids/fluorescent_protein_genes_and_plasmids/EGFP).
```python
egfp_sequence = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
with open("/scratch/s/shreejoy/nxu/Genomic_references/hg38/Raw/Homo_sapiens.GRCh38.cdna.ncrna.fa", "+a") as f:
    f.write(f">eGFP\n{egfp_sequence}")
```
3. From Genomic_references directory, generate salmon index
```bash
salmon --no-version-check index -t Raw/Homo_sapiens.GRCh38.cdna.ncrna.fa -i salmon_index_trap_seq
```
4. From project directory (`/scratch/s/shreejoy/nxu/trap_seq`), generate salmon scripts.  These commands have to be ran from the project directory since the scripts contain relative paths.
```bash
scripts/generate_salmon_scripts.sh data/WEN26595.20240222/240214_A02072_0139_AH2M75DSXC salmon /scratch/s/shreejoy/nxu/Genomic_references/hg38/salmon_index_trap_seq
```
5. Submit each salmon sbatch scripts to Niagara.
```
find salmon/salmon_scripts/ -type f | xargs -I {} sbatch {}
```