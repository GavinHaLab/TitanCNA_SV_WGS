# Snakemake workflow for TitanCNA_SV_WGS
#

# Modules to load
 * ml snakemake/5.19.2-foss-2019b-Python-3.7.4
 * ml R/3.6.2-foss-2019b-fh1
 * ml Python/3.7.4-foss-2019b-fh1
 * ml BCFtools/1.9-GCC-8.3.0
 * ml Pysam/0.15.4-GCC-8.3.0-Python-3.7.4
 * ml PyYAML/5.1.2-GCCcore-8.3.0-Python-3.7.4

# Set-up
## config/samples.yaml
Please specify the samples to be analyzed in config/samples.yaml, following the format explained therein.
 
## config/config.yaml
There are a number of parameters to adjust in config/config.yaml.  Filepaths to where your TitanCNA and ichorCNA repository as well as the filepath to tools (samTools, bcfTools, svaba) and readCounterScript.

* TitanCNA repository : https://github.com/gavinha/TitanCNA

* ichorCNA repository : https://github.com/GavinHaLab/ichorCNA

# Running the snakemake workflows on slurm cluster

This workflow will run TitanCNA (CN) anlaysis for a set of tumor-normal pairs, starting from the aligned BAM files. 

`snakemake -s TitanCNA.snakefile --latency-wait 60 --restart-times 3 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 30`

This workflow will run SvABA structural variation (SV) analysis for a set of tumor-normal pairs, starting from the aligned BAM files. 

`snakemake -s svaba.snakefile --latency-wait 60 --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 30`

SV classes are predicted by combining SV results with copy number results previously generated by TitanCNA to determine SV classes.

`snakemake -s combineSvabaTitan.snakefile --latency-wait 60 --keep-going  --restart-times 3 --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 30`
