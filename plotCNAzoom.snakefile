configfile: "config/configPlotZoom.yaml"
configfile: "config/samples_all.yaml"

'''
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml R/3.6.2-foss-2019b-fh1
ml Python/3.7.4-foss-2019b-fh1
ml BCFtools/1.9-GCC-8.3.0
ml Pysam/0.15.4-GCC-8.3.0-Python-3.7.4
ml PyYAML/5.1.2-GCCcore-8.3.0-Python-3.7.4

snakemake -s plotCNAzoom.snakefile --latency-wait 10 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100
'''

import glob
import os

def getTITANpath(id, ext):
    pattern = os.path.join(config["titanPath"], f"{id}_cluster*{ext}")
    matches = glob.glob(pattern)
    if not matches:
        raise ValueError(
            f"No TITAN file found for sample '{id}' with extension '{ext}'.\n"
            f"Searched: {pattern}"
        )
    return matches[0]


rule all:
  input: 
    expand("results/plotTITAN_zoom/{plotID}/{tumor}_CNA_{type}_chr{chr}-{start}-{end}.{format}",
           tumor=config["pairings"],
           plotID=config["plot_id"],
           type=config["plot_type"],
           chr=config["plot_chr"],
           start=config["plot_startPos"],
           end=config["plot_endPos"],
           format=config["plot_format"])


rule plotTITAN:
  input:
    titanBinFile=lambda wildcards: getTITANpath(wildcards.tumor, ".titan.ichor.cna.txt"),
    titanSegFile=lambda wildcards: getTITANpath(wildcards.tumor, ".titan.ichor.seg.noSNPs.txt"),
    titanParamFile=lambda wildcards: getTITANpath(wildcards.tumor, ".params.txt")
  output:
    "results/plotTITAN_zoom/{plotID}/{tumor}_CNA_{type}_chr{chr}-{start}-{end}.{format}"
  params:
    plotCNAscript=config["plotSVCN_script"],
    plotfuncs=config["plot_funcs"],
    libdir=config["titan_libdir"],
    genomeBuild=config["genomeBuild"],
    genomeStyle=config["genomeStyle"],
    cytobandFile=config["cytobandFile"],
    zoom=config["plot_zoom"],
    chr=config["plot_chr"],
    start=config["plot_startPos"],
    end=config["plot_endPos"],
    ylim=config["plot_ylim"],
    type=config["plot_type"],
    geneFile=config["plot_geneFile"],
    size=config["plot_size"],
    format=config["plot_format"]
  log:
    "logs/plotTITAN_zoom/{plotID}/{tumor}_{type}_chr{chr}-{start}-{end}.{format}.log"
  shell:
    "Rscript {params.plotCNAscript} "
    "--id {wildcards.tumor} "
    "--plot_funcs {params.plotfuncs} "
    "--titan_libdir {params.libdir} "
    "--titanBinFile {input.titanBinFile} "
    "--titanSegFile {input.titanSegFile} "
    "--titanParamFile {input.titanParamFile} "
    "--chrs \"{params.chr}\" "
    "--genomeBuild {params.genomeBuild} "
    "--genomeStyle {params.genomeStyle} "
    "--cytobandFile {params.cytobandFile} "
    "--start {params.start} "
    "--end {params.end} "
    "--zoom {params.zoom} "
    "--plotYlim \"{params.ylim}\" "
    "--geneFile {params.geneFile} "
    "--plotCNAtype {params.type} "
    "--plotSize \"{params.size}\" "
    "--outPlotFile {output} "
    "> {log} 2>&1"
