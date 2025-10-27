"""
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml R/3.6.2-foss-2019b-fh1
ml Python/3.7.4-foss-2019b-fh1
ml BCFtools/1.9-GCC-8.3.0
ml Pysam/0.15.4-GCC-8.3.0-Python-3.7.4
ml PyYAML/5.1.2-GCCcore-8.3.0-Python-3.7.4
ml annovar


#command to run snakemake (remove -np at end when done validating):
snakemake -s svaba_annovar.snakefile --latency-wait 60 --restart-times 2 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"
"""

configfile: "config/config_hg38.yaml"
#configfile: "config/samples_svaba_anno.yaml"
configfile: "config/samples_left.yaml"

rule all:
  input:
    expand("results/svaba/{tumor}/{tumor}.svaba.somatic.indel.multianno.vcf", tumor=config["pairings"])


rule run_annovar_indels:
	input:
		input_vcfs = "results/svaba/{tumor}/{tumor}.svaba.somatic.indel.vcf"
	output:
		output_vcf = "results/svaba/{tumor}/{tumor}.svaba.somatic.indel.multianno.vcf"
	params:
		annovar_python_script = config["annovar_python_script"],
		interpreter=config["interpreter"]
	log:
		"logs/runAnnovar_svaba/{tumor}_svaba_indels_runAnnovar.txt"
	shell:
		"({params.interpreter} {params.annovar_python_script} --input_vcf_file_path {input.input_vcfs}) 2> {log}"
