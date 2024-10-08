"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml R/3.6.2-foss-2019b-fh1
ml Python/3.7.4-foss-2019b-fh1
ml BCFtools/1.9-GCC-8.3.0
ml Pysam/0.15.4-GCC-8.3.0-Python-3.7.4
ml PyYAML/5.1.2-GCCcore-8.3.0-Python-3.7.4
ml ImageMagick/7.1.0-53-GCCcore-12.2.0

#command to run snakemake (remove -np at end when done validating):
snakemake -s TitanCNA.snakefile --latency-wait 60 --restart-times 3 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 20 -np
"""

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

include: "ichorCNA.snakefile"
include: "getAlleleCounts.snakefile"
import os.path

CLUST = {1:[1], 2:[1,2], 3:[1,2,3], 4:[1,2,3,4], 5:[1,2,3,4,5], 6:[1,2,3,4,5,6], 7:[1,2,3,4,5,6,7], 8:[1,2,3,4,5,6,7,8], 9:[1,2,3,4,5,6,7,8,9], 10:[1,2,3,4,5,6,7,8,9,10]}
PLOIDY = {2:[2], 3:[2,3], 4:[2,3,4]}


rule all:
	input: 
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		"results/titan/hmm/optimalClusterSolution.txt",
		"results/titan/hmm/optimalClusterSolution/",
		expand("results/CurationFiles/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.pdf", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]])

		
rule runTitanCNA:
	priority: 10
	input:
		alleleCounts="results/titan/tumCounts/{tumor}.tumCounts.txt",
		corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt"		
	output:		
		titan="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		param="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		segTxt="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
		seg="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
	params:
		outRoot="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
		titanRscript=config["TitanCNA_rscript"],
		libdir=config["TitanCNA_libdir"],
		numCores=config["TitanCNA_numCores"],
		normal=config["TitanCNA_normalInit"],
		chrs=config["TitanCNA_chrs"],
		sex=config["sex"],
		genomeStyle=config["genomeStyle"],
		genomeBuild=config["genomeBuild"],
		cytobandFile=config["cytobandFile"],
		estimatePloidy=config["TitanCNA_estimatePloidy"],
		estimateClonality=config["TitanCNA_estimateClonality"],
		estimateNormal=config["TitanCNA_estimateNormal"],
		centromere=config["centromere"],
		alphaK=config["TitanCNA_alphaK"],
		#alphaR=config["TitanCNA_alphaR"],
		#alleleModel=config["TitanCNA_alleleModel"],
		txnExpLen=config["TitanCNA_txnExpLen"],
		plotYlim=config["TitanCNA_plotYlim"]
	log:
		"logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.log"
	shell:
		"Rscript {params.titanRscript} --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {params.outRoot} --libdir {params.libdir} --id {wildcards.tumor} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --cytobandFile {params.cytobandFile} --chrs \"{params.chrs}\" --gender {params.sex} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --alphaK {params.alphaK} --txnExpLen {params.txnExpLen} --plotYlim \"{params.plotYlim}\" > {log} 2> {log}"
	
rule combineTitanAndIchorCNA:
	priority: 3
	input:
		titanSeg="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt", 
		titanBin="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		titanParam="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		ichorSeg="results/ichorCNA/{tumor}/{tumor}.seg.txt",
		ichorBin="results/ichorCNA/{tumor}/{tumor}.cna.seg",
		ichorParam="results/ichorCNA/{tumor}/{tumor}.params.txt"
	output:
		segFile="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt",
		binFile="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt",
	params:
		combineScript=config["TitanCNA_combineTitanIchorCNA"],
		libdir=config["TitanCNA_libdir"],
		centromere=config["centromere"],
		sex=config["sex"],
		mergeIchorHOMD=config["mergeIchorHOMD"],
		isPDXorCellLine=config["isPDXorCellLine"],
		correctSegmentsInBins=config["correctSegmentsInBins"]
	log:
		"logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.combineTitanIchorCNA.log"
	shell:
		"Rscript {params.combineScript} --libdir {params.libdir} --titanSeg {input.titanSeg} --titanBin {input.titanBin} --titanParam {input.titanParam} --ichorSeg {input.ichorSeg} --ichorBin {input.ichorBin} --ichorParam {input.ichorParam} --mergeIchorHOMD {params.mergeIchorHOMD} --isPDXorCellLine {params.isPDXorCellLine} --correctSegmentsInBins {params.correctSegmentsInBins} --sex {params.sex} --outSegFile {output.segFile} --outBinFile {output.binFile} --centromere {params.centromere} > {log} 2> {log}"	
	
rule selectSolution:
	priority: 2
	input:
		#ploidyDirs=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		resultFiles=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]])
	output:
		"results/titan/hmm/optimalClusterSolution.txt"
	params:
		solutionRscript=config["TitanCNA_selectSolutionRscript"],
		threshold=config["TitanCNA_solutionThreshold"]
	log:
		"logs/titan/selectSolution.log"
	shell:
		"""
		ploidyRun2=results/titan/hmm/titanCNA_ploidy2/
		if [ -d results/titan/hmm/titanCNA_ploidy3/ ]; then
			ploidyRun3=results/titan/hmm/titanCNA_ploidy3/
		else
			ploidyRun3=NULL
		fi
		if [ -d results/titan/hmm/titanCNA_ploidy4/ ]; then
			ploidyRun4=results/titan/hmm/titanCNA_ploidy4/
		else
			ploidyRun4=NULL
		fi
		Rscript {params.solutionRscript} --ploidyRun2 $ploidyRun2 --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output} > {log} 2> {log}
		"""
		
rule copyOptSolution:
	priority: 1
	input:
		"results/titan/hmm/optimalClusterSolution.txt"
	output:
		directory("results/titan/hmm/optimalClusterSolution/")
	params:
	log:
		"logs/titan/hmm/optSolution/copyOptSolution.log"
	shell:
		"""
		curDir=`pwd`
		for i in `cut -f11 {input} | grep -v "path"`;
		do
			echo -e "Copying $curDir/${{i}} to {output}"
			cp -r ${{curDir}}/${{i}}* {output}
		done		
		"""

rule createCurationFiles:
	priority: -1
	input:
		params_file="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
	output:
		"results/CurationFiles/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.pdf"
	params:
		script=config["TitanCNA_curation_script"]
	log:
		"logs/titanCuration/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}_curation.log"
	shell:
		"(python3 {params.script} --input_params_file_path {input.params_file} ) 2> {log}"


ruleorder: runTitanCNA > combineTitanAndIchorCNA > selectSolution > copyOptSolution > createCurationFiles
