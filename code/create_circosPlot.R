#' create_circosPlot.R
#' author: Minjeong Ko
#' Fred Hutchinson Cancer Research Center
#' contact: <mko@fredhutch.org>
#' date:  May 17, 2022
#' description: Generate circos plots of copy number and SV breakpoint arcs.

library(RCircos)
library(data.table)
library(optparse)
#data(RCircos.Gene.Label.Data);


option_list <- list(
  make_option(c("--id"), type = "character", help = "Sample ID"),
  make_option(c("--svFile"), type = "character", help = "Combined Svaba and Titan SV file (*svabaTitan.sv.PONfilter.bedpe)."),
  make_option(c("--cnFile"), type = "character", help = "Combined Svaba and Titan CNA file (*svabaTitan.cn.txt)."),  
  make_option(c("--geneFile"), type="character", help="Gene file for annotation (tab delimited, Chromosome, chromStart, chromEnd, Gene columns are required)"),
  make_option(c("--genomeBuild"), type="character", default="hg38", help="Geome build. Default: [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "UCSC", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--excludeCNoverlapType"), type = "character", default = "c(\"Unknown-ShortSVwithCN\")", help = "Exclude SV class."),
  make_option(c("--excludeTools"), type = "character", default = "c()", help = "Exclude SVs predicted by these tools"),
  make_option(c("--outPlotFile"), type="character", help="Path to output figure file.")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

id <- opt$id
svFile <- opt$svFile
cnFile <- opt$cnFile
geneFile <- opt$geneFile
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
excludeCNoverlapType <- eval(parse(text = opt$excludeCNoverlapType))
excludeTools <- eval(parse(text = opt$excludeTools))
outPlotFile <- opt$outPlotFile

if (!is.null(geneFile)) gene.lable.data <- read.table (geneFile, sep="\t", header=T)

#Load chromosome cytoband data
if(genomeBuild =="hg38"){
	data(UCSC.HG38.Human.CytoBandIdeogram)
	cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
}else{
	data(UCSC.HG19.Human.CytoBandIdeogram)
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
}
cyto.info <- as.data.table(cyto.info)

#Setup RCircos core components
RCircos.Set.Core.Components(cyto.info, tracks.inside=4, tracks.outside=0) 

# Can update params here
# params <- RCircos.Get.Plot.Parameters() #$char.width
# params$text.size <- 0.4
# params$track.padding  <- 0.07
# RCircos.Reset.Plot.Parameters(params)

# load SV data and filter
sv.data <- fread(svFile, head=T)
if (!is.null(excludeCNoverlapType)){
	sv.data <- sv.data[!CN_overlap_type %in% excludeCNoverlapType]
}
if (!is.null(excludeTools)){
	sv.data <- sv.data[!Tool %in% excludeTools]
}

#check LongRanger coordinates within cyto.info
#if coord is larger than cyto.info, then set it to the largest in cyto.info
chrEnds <- cyto.info[, max(chromEnd), by = Chromosome]
for (i in sv.data[, unique(chrom1)]){
	chrEndCoord <- chrEnds[Chromosome == i, V1]
	sv.data[chrom1 == i & start1 > chrEndCoord, start1 := chrEndCoord]
	sv.data[chrom1 == i & end1 > chrEndCoord, end1 := chrEndCoord]
	sv.data[chrom2 == i & start2 > chrEndCoord, start2 := chrEndCoord]
	sv.data[chrom2 == i & end2 > chrEndCoord, end2 := chrEndCoord]
}

# load CNV data
cnv.data <- fread(cnFile, head=T)
cnv.data.trim <- cnv.data[, c("Chromosome","Start","End","Corrected_Copy_Number")]


if (!is.null(geneFile)){
	pdf(outPlotFile, onefile=T, height=8, width=8, compress=TRUE)
	RCircos.Set.Plot.Area()
	RCircos.Chromosome.Ideogram.Plot()
	RCircos.Gene.Connector.Plot(gene.lable.data, track.num = 1, side = "in")
	RCircos.Gene.Name.Plot(gene.lable.data, name.col = 4, track.num = 2)
	RCircos.Histogram.Plot(cnv.data.trim, data.col=4, track.num=3, side="in")
	RCircos.Link.Plot(sv.data, track.num=4, by.chromosome=FALSE)
	dev.off()
} else{
	pdf(outPlotFile, onefile=T, height=8, width=8, compress=TRUE)
	RCircos.Set.Plot.Area()
	RCircos.Chromosome.Ideogram.Plot()
	RCircos.Histogram.Plot(cnv.data.trim, data.col=4, track.num=1, side="in")
	RCircos.Link.Plot(sv.data, track.num=2, by.chromosome=FALSE)
	dev.off()
}
