#' combineTITAN-ichor.R
#' author: Gavin Ha 
#' institution: Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date: July 23, 2018

#' requires R-3.3+
#' @import data.table
#' @import GenomicRanges
#' @import stringr
#' @import optparse

library(optparse)
library(TitanCNA)
library(stringr)
library(data.table)
library(GenomicRanges)
library(dplyr)

option_list <- list(
  make_option(c("--titanSeg"), type="character", help="TitanCNA segs.txt file. Required."),
	make_option(c("--titanBin"), type="character", help="TitanCNA titan.txt file. Required."),
	make_option(c("--titanParams"), type="character", help="TitanCNA params.txt file. Required."),
	make_option(c("--ichorSeg"), type="character", help="ichorCNA segs.txt file. Required."),
	make_option(c("--ichorBin"), type="character", help="ichorCNA cna.seg file. Required."),
	make_option(c("--ichorParams"), type="character", help="ichorCNA params.txt file. Required."),
	make_option(c("--ichorNormPanel"), type="character", default=NULL, help="Panel of normals; bin-level black list."),
	make_option(c("--mergeIchorHOMD"), type="logical", default = FALSE, help="Merge ichorCNA HOMD segment into final combined segments."),
	make_option(c("--isPDXorCellLine"), type="logical", default = FALSE, help="Insert TITAN missing hetSite data using corresponding bins from ichorCNA for specific samples"),
	make_option(c("--correctSegmentsInBins"), type="logical", default = FALSE, help="Correct segment copy number for consecutive bins that match in copy number but differ from segment copy number. Use only with --isPDXorCellLine"),
	make_option(c("--sex"), type="character", default="female", help="female or male. Default [%default]."),
	make_option(c("--libdir"), type="character", help="TitanCNA directory path to source R files if custom changes made."),
	make_option(c("--codedir"), type="character", help="TitanCNA coding directory path for downstream processes."),
	make_option(c("--outSegFile"), type="character", help="New combined segment file. Required"),
	make_option(c("--outBinFile"), type="character", help="New combined bin-level file. Required"),
	make_option(c("--centromere"), type="character", default=NULL, help="Centromere table.")
  )

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

titanSeg <- opt$titanSeg
titanBin <- opt$titanBin
titanParams <- opt$titanParams
ichorSeg <- opt$ichorSeg
ichorBin <- opt$ichorBin
ichorParams <- opt$ichorParams
ichorNormPanel <- opt$ichorNormPanel
mergeIchorHOMD <- as.logical(opt$mergeIchorHOMD)
rescueTITAN <- as.logical(opt$isPDXorCellLine)
correctSegmentsInBins <- as.logical(opt$correctSegmentsInBins)
gender <- opt$sex
outSegFile <- opt$outSegFile
outBinFile <- opt$outBinFile
centromere <- opt$centromere
libdir <- opt$libdir
codedir <- opt$codedir
outImageFile <- gsub(".seg.txt", ".RData", outSegFile)

# Get the output path for the repaired PDF file, using the same directory as the Bin file output path
if (rescueTITAN){
	outplot_repaired <- paste0(gsub("/[^/]*$", "", outBinFile), "/")
}

if (!is.null(libdir) && libdir != "None"){
  source(paste0(libdir, "/R/utils.R"))
  source(paste0(libdir, "/R/plotting.R"))
}

options(stringsAsFactors=F, width=150, scipen=999)
save.image(outImageFile)
## copy number state mappings ##
ichorCNmap <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP", rep("HLAMP", 1000))
#ichorCNmap <- list("0"="HOMD", "1"="DLOH", "2"="NEUT", "3"="GAIN", "4"="AMP", "5"="AMP")
maxichorcn <- 5

## load segments 
titan <- fread(titanSeg)
if (length(titan[["Cellular_Frequency"]] == 0)){
	setnames(titan, "Cellular_Frequency", "Cellular_Prevalence")
}
ichor.segs <- fread(ichorSeg)
setnames(ichor.segs, c("ID", "chrom", "start", "end", "num.mark", "seg.median.logR", "copy.number", "call"), 
		c("Sample", "Chromosome", "Start_Position.bp.", "End_Position.bp.", 
		  "Length.snp.", "Median_logR", "Copy_Number", "TITAN_call"))

## load data points ##
titan.cn <- fread(titanBin)
titan.cn <- cbind(Sample=titan[1,Sample], titan.cn)
#titan.cn[, chr := as.character(Chr)]
id <- titan[1, Sample]
ichor.cn <- fread(ichorBin)
ichor.cn <- cbind(Sample = id, ichor.cn)
#ichor.cn[, CopyNumber := state - 1]
ichor.cn[, Position := start]
setnames(ichor.cn, c("chr", "start", paste0(id,".copy.number"), paste0(id,".event"), paste0(id,".logR"), "end"), 
		c("Chr", "Start", "CopyNumber", "TITANcall", "LogRatio", "End"))

## get chromosome style
titan$Chromosome <- as.character(titan$Chromosome)
titan.cn$Chr <- as.character(titan.cn$Chr)
genomeStyle <- seqlevelsStyle(titan.cn$Chr)[1]
chrs <- c(1:22, "X")
chrXStr <- grep("X", chrs, value=TRUE)
seqlevelsStyle(chrs) <- genomeStyle

## load parameters ##
params <- read.delim(titanParams, header=F, as.is=T)
purity <- 1 - as.numeric(params[1,2])
ploidyT <- as.numeric(params[2,2])
ploidy <- purity * ploidyT + (1-purity) * 2
params.ichor <- read.delim(ichorParams, header=T, as.is=T)
homd.var <- as.numeric(strsplit(params[params[,1]=="logRatio Gaussian variance:",2], " ")[[1]][1])
homd.sd <- sqrt(homd.var)

## get gender
if (is.null(gender) || gender == "None"){
	gender <- params.ichor[3, 2]
}

## get bin overlap with SNPs - include ichor bins even if no SNPs overlap 
titan.gr <- titan.cn[, .(Chr, Position)]
titan.gr[, Start := Position]; titan.gr[, End := Position]
titan.gr <- as(titan.gr, "GRanges")
ichor.gr <- as(ichor.cn, "GRanges")
hits <- findOverlaps(query = titan.gr, subject = ichor.gr)
titan.cn[queryHits(hits), Start := ichor.cn[subjectHits(hits), Start]]
titan.cn[queryHits(hits), End := ichor.cn[subjectHits(hits), End]]
titan.ichor.cn <- merge(titan.cn, ichor.cn, by=c("Sample", "Chr", "Start", "End"), all=T, suffix=c("",".ichor"))
titan.ichor.cn[is.na(LogRatio), LogRatio := LogRatio.ichor] # assign ichor log ratio to missing titan SNPs
titan.ichor.cn <- titan.ichor.cn[, -c(grep("ichor", colnames(titan.ichor.cn),value=T)), with=F] 

## include HOMD from ichorCNA if missing in TitanCNA (usually when no SNPs are there)
if (mergeIchorHOMD){
	homdLogRThres.auto <- log2((2*(1-purity)) / (2*(1-purity) + ploidyT*purity)) + homd.sd
	ichor.segs[, logR_Copy_Number := logRbasedCN(Median_logR, purity, ploidy, cn=2)]
	ichor.segs[, Corrected_Copy_Number := as.integer(round(logR_Copy_Number))]
	ichor.segs[, Corrected_Call := ichorCNmap[Corrected_Copy_Number + 1]]
	ichor.segs[, Corrected_logR := log2(logR_Copy_Number / ploidy)]
	ichor.segs.homd.ind <- ichor.segs[Chromosome != chrXStr & Corrected_Copy_Number == 0 & Median_logR < homdLogRThres.auto, which = TRUE]
	if (length(ichor.segs.homd.ind)){
		titan.gr <- copy(titan)
		titan.gr[, Start := Start_Position.bp.]; titan.gr[, End := End_Position.bp.]
		titan.gr <- as(titan.gr, "GRanges")
		ichor.segs.gr <- copy(ichor.segs)
		ichor.segs.gr[, Start := Start_Position.bp.]; ichor.segs.gr[, End := End_Position.bp.]
		ichor.homd.gr <- as(ichor.segs.gr[ichor.segs.homd.ind], "GRanges")
		#hits <- findOverlaps(query = titan.gr, subject = ichor.homd.gr)
		# chromosome X automatically ignored since titan doesn't contain this seqname
		titan.gr.combinedHomd <- sort(c(titan.gr, ichor.homd.gr))
		titan.combinedHomd <- as.data.table(disjoin(titan.gr.combinedHomd, with.revmap = TRUE))
		mcolData <- data.table()
		for (i in 1:nrow(titan.combinedHomd)){
			ind <- titan.combinedHomd[i, revmap][[1]]
			if (length(ind) == 1){ # same titan segment with no split
				mcolData <- rbind(mcolData, mcols(titan.gr.combinedHomd)[ind, c(4:19,21)])
			}else{ # ind has more than one mapping index
				ind.homd.gr <- ind[which(titan.gr.combinedHomd[ind]$Corrected_Copy_Number == 0)] 
				# should only be single index value
				if (length(ind.homd.gr) > 1){ # both indices are HOMD, then choose titan results
					ind.homd.gr <- ind[which(!is.na(titan.gr.combinedHomd$Median_Ratio[ind]))]
				}
				mcolData <- rbind(mcolData, mcols(titan.gr.combinedHomd)[ind.homd.gr[1], c(4:19,21)])
			}
		}
		titan.combinedHomd <- cbind(titan.combinedHomd, as.data.table(mcolData))
		titan <- cbind(Sample = id, titan.combinedHomd[, -c(4:6)])
		setnames(titan, c("seqnames", "start", "end"), c("Chromosome", "Start_Position.bp.", "End_Position.bp."))
	}
}

## combine TITAN (chr1-22) and ichorCNA (chrX) segments and bin/SNP level data ##
## if male only ##
if (gender == "male"){
	cn <- rbind(titan.ichor.cn[Chr %in% chrs[1:22]], ichor.cn[Chr == chrs[grep("X", chrs)]], fill = TRUE)
	segs <- rbind(titan[Chromosome %in% chrs[1:22]], ichor.segs[Chromosome == chrs[grep("X", chrs)]], fill = TRUE)
	segs[, subclone.status := NULL]
}else{
	cn <- titan.ichor.cn
	segs <- titan
}

## sort column order
setnames(segs, c("Start_Position.bp.", "End_Position.bp."), c("Start", "End"))
cols <- c("Sample", "Chr", "Position", "Start", "End")
setcolorder(cn, c(cols, colnames(cn)[!colnames(cn) %in% cols]))

## get major/minor CN from segs and place in SNP/level data ##
#cn.gr <- cn[, .(Chr, Start, End)]
#cn.gr <- as(na.omit(cn.gr), "GRanges")
#segs.gr <- as(segs, "GRanges")
#hits <- findOverlaps(query = cn.gr, subject = segs.gr)
#cn[queryHits(hits), MajorCN := segs[subjectHits(hits), MajorCN]]
#cn[queryHits(hits), MinorCN := segs[subjectHits(hits), MinorCN]]
#cn[is.na(CopyNumber), CopyNumber := MajorCN + MinorCN]

## remove LogRatio for bins that do not overlap segments (i.e. remaining rows with CopyNumber == NA
# these regions that are not part of segments are usually noisy
#cn[is.na(CopyNumber), LogRatio := NA]

## filter outlier HOMD data points
# message("Filtering bins with outlier negative log ratios...")
# homdLenThres <- 10000
# homdNumSNPThres <- 40
# homdLogRThres.auto <- log2((2*(1-purity)) / (2*(1-purity) + ploidyT*purity)) - 1*homd.sd
# ## OUTLIERS IN CHROMOSOME X SHOULD BE AGNOISTIC OF SEX 
# #if (gender == "male"){
# #	homdLogRThres.X <- round(log2((1*(1-purity)) / (1*(1-purity) + ploidyT*purity)), digits = 1) - 0.1
# #}else{
# #	homdLogRThres.X <- homdLogRThres.auto
# #}
# ind.homd.cn <- cn[LogRatio < homdLogRThres.auto, which = TRUE]
# message("Removing ", length(ind.homd.cn), " negative log ratio outlier bins ...")
# ind.segs.remove <- segs[((End-Start < homdLenThres | Length.snp. < homdNumSNPThres) & Corrected_Call == "HOMD") | Median_logR < homdLogRThres.auto, which=T]
# message("Removing ", length(ind.segs.remove), " negative log ratio outlier segments ...")
# if (length(ind.segs.remove) > 0){
# 	segsToRemove <- segs[ind.segs.remove]
# 	hits <- findOverlaps(query = as(as.data.frame(segsToRemove), "GRanges"), subject = as(as.data.frame(cn), "GRanges"))
# 	ind.homd.segs <- union(subjectHits(hits), ind.homd.cn)
# }else{
# 	ind.homd.segs <- ind.homd.cn
# }

#message("Loading panel of normals: ", ichorNormPanel)
#panel <- readRDS(ichorNormPanel)
#panel.dt <- as.data.table(as.data.frame(panel))
#panel.colNames <- colnames(panel.dt[, 11:(ncol(panel.dt)-1)])
#panel.dt[, variance := apply(.SD, 1, var, na.rm=TRUE), .SDcols = panel.colNames]
#panel.sd <- 2 * sd(panel.dt$variance, na.rm=T)
#panel.dt[, outlier := variance < -panel.sd | variance > panel.sd]
#ind.panel <- which(panel$Median < -2*sd(x$panel,na.rm=T))
#ind.panel <- which(panel.dt$outlier == TRUE)
#hits <- findOverlaps(query = as(as.data.frame(cn), "GRanges"), subject = as(panel[ind.panel,], "GRanges"))
#ind.homd.panel <- queryHits(hits)
#ind.homd.all <- union(ind.homd.panel, ind.homd.segs)
# remove the elements 
# if (length(ind.homd.segs) > 0){
# 	cn <- cn[-ind.homd.segs]
# }
# #cn <- cn[ind.homd.segs, FILTER := "EXCLUDE"]
# if (length(ind.segs.remove) > 0){
# 	segs <- segs[-ind.segs.remove]
# }
#segs <- segs[ind.segs.remove, FILTER := "EXCLUDE"]

if (rescueTITAN) {
	logRCN_name <- paste0(id,".logR_Copy_Number")
	CCN_name <- paste0(id,".Corrected_Copy_Number")
	cn_repaired <- fill_cellular_prevalence(cn, logRCN_name, CCN_name)

	## correct copy number beyond maximum CN state based on purity and logR
	correctCN <- correctIntegerCN(cn_repaired, segs, purity, ploidyT, maxCNtoCorrect.autosomes = NULL, 
								maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.05, gender = gender, chrs = chrs)
	segs <- correctCN$segs
	cn <- correctCN$cn
} else {
	## correct copy number beyond maximum CN state based on purity and logR
	correctCN <- correctIntegerCN(cn, segs, purity, ploidyT, maxCNtoCorrect.autosomes = NULL, 
								maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.05, gender = gender, chrs = chrs)
	segs <- correctCN$segs
	cn <- correctCN$cn
}

centromeres <- fread(centromere)

#### For cases where sample is high-purity tumor ran in tumor-only mode ####
if (rescueTITAN == TRUE) {
  ##################################### FIXING SEGS FILE OUTPUT #####################################
  ## Adjust ichorCNA segment copy numbers to match TitanCNA baseline
  ichor_2.segs <- copy(ichor.segs)
  ichor_2.segs[grepl("X", Chromosome), logR_Copy_Number := logRbasedCN(Median_logR, purity, ploidyT, cn=2)]
  ichor_2.segs[!grepl("X", Chromosome), logR_Copy_Number := logRbasedCN(Median_logR, purity, ploidyT, cn=1)]
  ichor_2.segs[, Corrected_Copy_Number := as.integer(round(logR_Copy_Number))]
  ichor_2.segs[, Corrected_Call := ichorCNmap[Corrected_Copy_Number + 1]]
  ichor_2.segs[, Corrected_logR := log2(logR_Copy_Number / ploidy)]
  ichor_2.segs[, Start := Start_Position.bp.]; ichor_2.segs[, End := End_Position.bp.]
  
  ## Load file
  titan.gr <- copy(segs)
  
  ## Remove "sparse" regions from TITAN
  print("...repairing sparse regions from TITAN seg file")
  # Remove the individual rows in which there is only 1 SNP per 1-bp length segment
  titan_subset <- subset(titan.gr, (titan.gr$End != (titan.gr$Start))) # PRIMARY ISSUE: needed to remove TitanCNA's single position calls to accomodate ichor's segments
  # Create a vector that calculates the ratio of the length of the segment to the number of SNPs in the segment
  titan_subset$Density.snp. <- ((titan_subset$Length.snp.) / (titan_subset$End - titan_subset$Start))
  # Get the 0.05 quantile threshold for the density of SNPs in the segments
  density_quantile_05 <- quantile(titan_subset$Density.snp., probs=c(0.05))
  # Create a vector that annotates for each row that, if the density of SNPs in the segment is greater than 2 standard deviations from the mean, then it is a "sparse" segment - annotate sparse with -1, otherwise 0
  titan_subset$Sparse <- ifelse(titan_subset$Density.snp. < (density_quantile_05), -1, 0)
  
  ## If the segment is sparse and the length.snp. is less than 100, then remove the segment
  titan_subset_filtered <- subset(titan_subset, (titan_subset$Sparse != -1) | (titan_subset$Length.snp. > 100))
  
  titan_subset_filtered[, Start := Start]; titan_subset_filtered[, End := End]
  
  ########## Fix bin and segment files ##########
  print("...Fixing Bin File for Empty Areas")
  cn2 <- fill_cn_blanks_if_ichor(cn)
  
  ## Initialize segments for regions where gaps exceed 10kb
  titan_subset_filtered <- initialize_new_segments(titan_subset_filtered, ichor_2.segs)
  
  
  print("...Beginning sparsity removal re-fill of values")
  ## Fill in the new segments with median values from corresponding ichorCNA bins
  titan_subset_filtered <- fill_new_segments(titan_subset_filtered, cn2)
  
  
  print("...processing segments")
  segs2 <- extendSegments(titan_subset_filtered, removeCentromeres = TRUE, centromeres = centromeres, extendToTelomeres = FALSE,
                          chrs = chrs, genomeStyle = genomeStyle)
  
  segs_plotting <- copy(segs2)
  setnames(segs_plotting, c("Start", "End"), c("Start_Position.bp.", "End_Position.bp."))
  
  ########## Process cn_plotting ##########
  cn_plotting <- cn2[ (which(!is.na(cn2$Position) & !is.na(cn2$logR_Copy_Number))),]
  cn_plotting <- cn_plotting[,c("Chr", "Position", "RefCount", "Depth", "AllelicRatio", "LogRatio", "CopyNumber", "TITANstate", "TITANcall", "ClonalCluster", "CellularPrevalence", "logR_Copy_Number", "Corrected_logR", "Corrected_Ratio", "Corrected_Copy_Number", "Corrected_Call")]
  cn_plotting$Chr <- factor(cn_plotting$Chr, levels = chrs)
  cn_plotting <- cn_plotting[order(cn_plotting$Chr), ]

  ########## Process segment-based copy number ##########
  ## Create a new column for segment-based copy number in the cn_plotting category
  seg_copy_numbers <- numeric(length = nrow(cn_plotting))
  print("...processing segment_based copy")
  for (i in 1:nrow(cn_plotting)) {
    # Check if POSITION falls between START and END in table 2
    matched_row <- which(cn_plotting$Chr[i] == segs_plotting$Chromosome & cn_plotting$Position[i] >= segs_plotting$Start_Position.bp. & cn_plotting$Position[i] <= segs_plotting$End_Position.bp.)
    # If a match is found, extract the copy number
    if (length(matched_row) > 0) {
      seg_copy_numbers[i] <- segs_plotting$Corrected_Copy_Number[matched_row]
    } else {
      seg_copy_numbers[i] <- cn_plotting$Corrected_Copy_Number[i]  # If no match found, store Corrected_Copy_Number
    }
    if (cn_plotting$Chr[i] == "chrX") {
      cn_plotting$Corrected_Copy_Number[i] <- seg_copy_numbers[i]
    }
  }
  ## Assign segment copy numbers to the cn_plotting table
  cn_plotting$Segment_Copy_Number <- seg_copy_numbers
  
  # Optional fine-tuning: if multiple bins have a shared copy number that also differs from Segment_Copy_Number, replace Segment_Copy_Number with the shared copy number
  # Threshold represents the minimum number of bins that have to share copy number to be considered
  if (correctSegmentsInBins) {
	corCN_col = "Corrected_Copy_Number"
	segCN_col = "Segment_Copy_Number"
	cn_plotting <- correct_segCN_bins(cn_plotting, CorrectedCN_Col=corCN_col, SegCN_Col=segCN_col, threshold=50)
  }
  
  ########## Plotting repaired Titan-Ichor Genome-Wide Plots ##########
  print("...loading plotting parameters")
  
  chrs_plot <- c(chrs, "chrX")
  
  ## Plotting Titan-Ichor Repaired File
  cn_plotting_dt <- as.data.table(cn_plotting)
  
  # Load  in TitanSvaba plotting library to load plotTitanIchorCNA function
  source(paste0(codedir, "/plotting.R"))
  
  ## Invoke plotting parameters and function
  ("...preparing to plot")
#   normCN <- 2
  yaxis <- "integer"
  ylim <- c(-2, 4)
#   ploidyS <- purity * ploidyT + (1-purity) * normCN
#   ploidyX <- purity * ploidyT/2 + (1 - purity) * normCN
  if (yaxis == "integer"){
    cn_plotting_dt[!grepl("X",Chr), logR_Copy_Number := (logRbasedCN(LogRatio, purity, ploidyT, cn=2))]
    cn_plotting_dt[grepl("X",Chr), logR_Copy_Number := (logRbasedCN(LogRatio, purity, ploidyT, cn=1))]
    cn_plotting_dt[, Corrected_Copy_Number := round(logR_Copy_Number)]
    
    # Correct Corrected_Copy_Number values without Cellular_Prevalence using logRbasedCN within both titan_ichor_merged and titan_ichor_new_method
    # Create a new column in titan_ichor_merged called "Adjusted_Copy_Number" that uses logRbasedCN to correct the Copy_Number values
    segs_plotting[!grepl("X", Chromosome), Adjusted_Copy_Number := round(logRbasedCN(Median_logR, purity, ploidyT, cn=2))]
    segs_plotting[grepl("X", Chromosome), Adjusted_Copy_Number := round(logRbasedCN(Median_logR, purity, ploidyT, cn=1))]
    
    adj_copy_numbers <- numeric(length = nrow(cn_plotting_dt))
    
    for (i in 1:nrow(cn_plotting_dt)) {
      # Check if POSITION falls between START and END in table 2
      matched_row <- which(cn_plotting_dt$Chr[i] == segs_plotting$Chromosome & cn_plotting_dt$Position[i] >= segs_plotting$Start_Position.bp. & cn_plotting_dt$Position[i] <= segs_plotting$End_Position.bp.)
      # If a match is found, extract the copy number
      if (length(matched_row) > 0) {
        adj_copy_numbers[i] <- segs_plotting$Adjusted_Copy_Number[matched_row]
      } else {
        adj_copy_numbers[i] <- cn_plotting_dt$Corrected_Copy_Number[i]  # If no match found, store NA
        # adj_copy_numbers[i] <- NA
      }
      if (cn_plotting_dt$Chr[i] == "chrX") {
        cn_plotting_dt$Corrected_Copy_Number[i] <- adj_copy_numbers[i]
      }
    }
    cn_plotting_dt$Adjusted_Segment_Copy_Number <- adj_copy_numbers

    ### Adjust segment copy numbers again to account for new calculation without CP + generate new "segments" based on CN copy number distribution
	if (correctSegmentsInBins) {
		cn_plotting_dt2 <- copy(cn_plotting_dt)

		for (i in 1:nrow(cn_plotting_dt2)) {
			if (cn_plotting_dt2$Adjusted_Segment_Copy_Number[i] != cn_plotting_dt2$Corrected_Copy_Number[i]) {
				counter <- 1
				for (j in (i+1):nrow(cn_plotting_dt2)) {
				if (cn_plotting_dt2$Adjusted_Segment_Copy_Number[j] != cn_plotting_dt2$Corrected_Copy_Number[j]) {
					counter <- counter + 1
				} else {
					break
				}
				}
				if (counter > 10) {
					for (j in i:(i+counter)) {
						cn_plotting_dt2$Adjusted_Segment_Copy_Number[j] <- cn_plotting_dt2$Corrected_Copy_Number[j]
					}
				}
			}
		}
	}
    colName <- "logR_Copy_Number"
  }else{
    cn_plotting_dt[!grepl("X",Chr), LogRatio := LogRatio + log2(ploidyS / 2)]
    cn_plotting_dt[grepl("X", Chr), LogRatio := LogRatio + log2(ploidyX / 2)]
    segs_plotting$LogRatio <- segs_plotting$Median_logR
    segs_plotting$LogRatio <- segs_plotting$LogRatio + log2(ploidyS / 2)
    colName <- "LogRatio"
  }
  
  print("Plotting...")
  outFile_pdf <- paste0(outplot_repaired, id, "_TitanIchor_Repaired_GenomeWide_logRCN.pdf")
  pdf(paste0(outFile_pdf),width=20,height=6)
  plotTitle <- paste0(id, "_TitanIchorCNA GenomeWide LogR-Copy-Number (Adjusted Segment CN Coloring)")
  plotTitanIchorCNA(as.data.frame(cn_plotting_dt), chr=chrs_plot, colName=colName, callColName="Adjusted_Segment_Copy_Number", cytoBand=FALSE, purity=purity, ploidyT=ploidyT, yaxis=yaxis,
                    cnCol=NULL, yrange=ylim, genomewide=TRUE, spacing=4, xaxt="n", cex=0.5, gene.cex=1.5, plot.title=plotTitle)
  dev.off()
  
  cn_out <- copy(cn2)
  segs_out <- copy(segs2)
  
  write.table(cn_plotting_dt, file = paste0(outplot_repaired, id, ".TitanIchor.repaired-seg.cna.plotting.txt"), col.names=T, row.names=F, quote=F, sep="\t")
  write.table(segs_plotting, file = paste0(outplot_repaired, id, ".TitanIchor.repaired-seg.segs.plotting.txt"), col.names=T, row.names=F, quote=F, sep="\t")
  
} else {
  segs_out <- extendSegments(segs, removeCentromeres = TRUE, centromeres = centromeres, extendToTelomeres = FALSE,
                             chrs = chrs, genomeStyle = genomeStyle)
  cn_out <- copy(cn)
}

## write segments to file ##
write.table(segs_out, file = outSegFile, col.names=T, row.names=F, quote=F, sep="\t")
write.table(cn_out, file = outBinFile, col.names=T, row.names=F, quote=F, sep="\t")
## write segments without germline SNPs
outSegNoSNPFile <- gsub(".txt", ".noSNPs.txt", outSegFile)
write.table(segs_out[, -c("Start.snp", "End.snp")], file = outSegNoSNPFile, col.names=T, row.names=F, quote=F, sep="\t")

save.image(outImageFile)

## write segments in IGV / GISTIC format ##
igv <- segs_out[, .(Sample, Chromosome, Start.snp, End.snp, Length.snp., logR_Copy_Number)]
igv[Chromosome %in% chrs[1:22], Corrected.logR := log2(logR_Copy_Number / 2)]
igv[Chromosome == chrs[grep("X", chrs)], Corrected.logR := log2(logR_Copy_Number / 1)]
igv[, logR_Copy_Number := NULL]
outIGVFile <- gsub("seg.txt", "segIGV.txt", outSegFile)
write.table(igv, file = outIGVFile, col.names=T, row.names=F, quote=F, sep="\t")



