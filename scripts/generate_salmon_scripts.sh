#!/bin/bash

datadir="$1" #directory of where all the data/fastqfiles are stored (absolute path)
outputdir="$2" #directory of where you want the output to be stored (absolute path)
indexdir="$3" #directory of salmon index (absolute path)
num_threads=80

#initial setup/prep
echo "Processing files from $datadir ."

if [ ! -d "${outputdir}/salmon_scripts" ]; then
  mkdir -p "${outputdir}/salmon_scripts"
fi

if [ ! -d "${outputdir}/salmon_results" ]; then
  mkdir -p "${outputdir}/salmon_results"
fi

if [ ! -d "${outputdir}/slurm_logs" ]; then
  mkdir -p "${outputdir}/slurm_logs"
fi

res_fold="${outputdir}/salmon_results/"

arr=("${datadir}"/*_L004_R1_001.fastq.gz)
files_of_interest=( "${arr[@]/$'_L004_R1_001.fastq.gz'}" )
for i in "${!files_of_interest[@]}"; do
    files_of_interest[i]=${files_of_interest[i]#"${datadir}/"}
done

tempi=0

for file in "${files_of_interest[@]}"; do
    r1=${file}_L004_R1_001.fastq.gz
    r2=${file}_L004_R2_001.fastq.gz
    
    echo "Processing sample ${file} using $datadir/$r1 and $datadir/$r2"
    
    #Writing the script to process each sample
    {
      echo '#!/bin/bash'
      echo '#SBATCH --time=01:00:00'
      echo '#SBATCH --nodes=1'
      echo '#SBATCH --cpus-per-task=80'
      echo "#SBATCH --job-name=${file}"
      echo "#SBATCH --output=${outputdir}/slurm_logs/${file}.out"
      echo "conda activate /scratch/s/shreejoy/nxu/CIHR/envs"
      echo "salmon --no-version-check quant -i \"$indexdir\" -l A -1 \"${datadir}/$r1\" -2 \"${datadir}/$r2\" -p $num_threads -o \"${res_fold}${file}_quant\""
    } > "${outputdir}/salmon_scripts/SalmonParaScript${tempi}.sh"

    chmod +x "$outputdir/salmon_scripts/SalmonParaScript${tempi}.sh"
    echo "A script has been generated for $file . Script has the name STARParaScript$tempi.sh and its results will start with $file ."
    ((tempi=tempi+1))
done
