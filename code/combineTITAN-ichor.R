#' combineTITAN-ichor.R
#' author: Gavin Ha 
#' institution: Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date: July 23, 2018

# cur_path <- /Volumes/fh/fast/ha_g/projects/NelsonLab/MedGenome_WGS_LNCaP/matrices_CNA_06192023/TITAN_outputs/NEW_combineTITAN-ichor.R

#' requires R-3.3+
#' @import data.table
#' @import GenomicRanges
#' @import stringr
#' @import optparse

library(optparse)

option_list <- list(
  make_option(c("--titanSeg"), type="character", help="TitanCNA segs.txt file. Required."),
  make_option(c("--titanBin"), type="character", help="TitanCNA titan.txt file. Required."),
  make_option(c("--titanParams"), type="character", help="TitanCNA params.txt file. Required."),
  make_option(c("--ichorSeg"), type="character", help="ichorCNA segs.txt file. Required."),
  make_option(c("--ichorBin"), type="character", help="ichorCNA cna.seg file. Required."),
  make_option(c("--ichorParams"), type="character", help="ichorCNA params.txt file. Required."),
  make_option(c("--ichorNormPanel"), type="character", help="Panel of normals; bin-level black list."),
  make_option(c("--mergeIchorHOMD"), type="logical", default = FALSE, help="Merge ichorCNA HOMD segment into final combined segments."),
  make_option(c("--rescueTITANgaps"), type="logical", default = FALSE, help="Fill in TITAN missing hetSite data using corresponding bins from ichorCNA"),
  make_option(c("--sex"), type="character", default="female", help="female or male. Default [%default]."),
  make_option(c("--libdir"), type="character", help="TitanCNA directory path to source R files if custom changes made."),
  make_option(c("--outSegFile"), type="character", help="New combined segment file. Required"),
  make_option(c("--outBinFile"), type="character", help="New combined bin-level file. Required"),
  make_option(c("--centromere"), type="character", default=NULL, help="Centromere table."),
  make_option(c("--titan_Rdata"), type="character", default=NULL, help="Output RData from TitanCNA run for Sample")
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
rescueTITAN <- as.logical(opt$rescueTITANgaps)
gender <- opt$sex
outSegFile <- opt$outSegFile
outBinFile <- opt$outBinFile
centromere <- opt$centromere
libdir <- opt$libdir
titan_rdata <- opt$titan_Rdata
# outplots <- opt$outplots
# Get the output path for the repaired PDF file, using the same directory as the Bin file output path
outplot_repaired <- paste0(gsub("/[^/]*$", "", outBinFile), "/")
outImageFile <- gsub(".seg.txt", ".RData", outSegFile)


library(devtools)
library(TitanCNA)
library(stringr)
library(data.table)
library(GenomicRanges)
library(dplyr)

if (!is.null(libdir) && libdir != "None"){
  source(paste0(libdir, "/R/utils.R"))
  source(paste0(libdir, "/R/plotting.R"))
}

options(stringsAsFactors=F, width=150, scipen=999)
save.image(outImageFile)
## copy number state mappings ##
ichorCNmap <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP", rep("HLAMP", 1000))
# ichorCNmap <- list("0"="HOMD", "1"="HETD", "2"="NEUT", "3"="GAIN", "4"="AMP", "5"="AMP")
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

## correct copy number beyond maximum CN state based on purity and logR
correctCN <- correctIntegerCN(cn, segs, purity, ploidyT, maxCNtoCorrect.autosomes = NULL, 
                              maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.05, gender = gender, chrs = chrs)
segs <- correctCN$segs
cn <- correctCN$cn

centromeres <- fread(centromere)

##  ----------------------------- ADDED ----------------------------- ##

adjustPositionsForInsertions <- function(table) {
  for (row in 2:(nrow(table) - 1)) {
    cur_row <- table[row,]
    prev_row <- table[(row - 1),]
    next_row <- table[(row+1),]
    if ((is.na(cur_row$Length.snp.)) & (all(sapply(list(prev_row$Chromosome, cur_row$Chromosome, next_row$Chromosome), function(x) x == cur_row$Chromosome)))) {
      if (((cur_row$Start) < (prev_row$End))) {
        table[row, "Start"] <- ((table[row-1, "End"]) + 1)
      }
      if ((cur_row$End) > (next_row$Start)) {
        table[row, "End"] <- ((table[row+1, "Start"]) - 1)
      }
      if ((cur_row$End) < (cur_row$Start)) {
        table[row, "End"] <- ((table[row, "Start"]) + 1)
      }
    }
  }
  return(table)
}


fill_cn_blanks_if_ichor <- function(table) {
  table %>%
    mutate(Position = if_else(
      is.na(Position) & !is.na(.[[ncol(.)]]), 
      .$Start, 
      Position
    )) %>%
    return()
}

library(dplyr)

## Use SubsetByOverlaps to reconcile missing data from TITAN using ichorCNA results ##
if (rescueTITAN) {
  ##################################### FIXING SEGS FILE OUTPUT #####################################
  ## Fix ichorCNA segs
  ichor_2.segs <- copy(ichor.segs)
  ichor_2.segs[, logR_Copy_Number := logRbasedCN(Median_logR, purity, ploidy, cn=2)]
  ichor_2.segs[, Corrected_Copy_Number := as.integer(round(logR_Copy_Number))]
  ichor_2.segs[, Corrected_Call := ichorCNmap[Corrected_Copy_Number + 1]]
  ichor_2.segs[, Corrected_logR := log2(logR_Copy_Number / ploidy)]
  
  ## Load file
  titan.gr <- copy(segs)
  
  # FIRST, remove all sparse regions from TITAN
  titan_subset.gr <- subset(titan.gr, ((titan.gr$End != (titan.gr$Start)) & (titan.gr$Length.snp. / (titan.gr$End - titan.gr$Start))  >= 0.0002)) # PRIMARY ISSUE: needed to remove TitanCNA's single position calls to accomodate ichor's segments
  titan_subset_filtered <- copy(titan_subset.gr)
  
  titan_subset.gr[, Start := Start]; titan_subset.gr[, End := End]
  titan_2.gr <- as(titan_subset.gr, "GRanges")
  
  ichor.segs.gr <- copy(ichor_2.segs)
  ichor.segs.gr[, Start := Start_Position.bp.]; ichor.segs.gr[, End := End_Position.bp.]
  ichor.segs.gr <- as(ichor.segs.gr, "GRanges")
  
  ## subsetByOverlaps - NOT Overlaps
  nohits <- subsetByOverlaps(ichor.segs.gr, titan_2.gr, type="any", minoverlap=10000, invert=TRUE)
  sbo_no_df <- DataFrame(nohits)
  
  sbo_no_df_list <- as.list(sbo_no_df)
  sbo_no_df_data <- as.data.frame(sbo_no_df_list$X)
  sbo_no_df_data <- subset(sbo_no_df_data, select=-c(start, end, width, strand))
  
  titan_ichor_merged <- merge(titan_subset.gr, sbo_no_df_data, by.x=c("Sample", "Chromosome", "Start", "End"), by.y=c("Sample", "seqnames", "Start_Position.bp.", "End_Position.bp."), all=T, suffix=c("",".ichor"))
  
  titan_ichor_merged <- adjustPositionsForInsertions(titan_ichor_merged)
  
   ##################################### Filling further GAPS  #####################################
  ## Fixing Copy Number Values
  titan_ichor_merged[is.na(Corrected_Call), Corrected_Call := Corrected_Call.ichor] # assign ichor Corrected_Calls to missing titan SNPs
  titan_ichor_merged[is.na(Corrected_Copy_Number), Corrected_Copy_Number := Corrected_Copy_Number.ichor] # assign ichor Corrected_copy_numbers to missing titan SNPs
  titan_ichor_merged[is.na(Corrected_MinorCN), Corrected_MinorCN := 0]
  titan_ichor_merged[is.na(Corrected_MajorCN), Corrected_MajorCN := Corrected_Copy_Number.ichor]
  titan_ichor_merged[is.na(logR_Copy_Number), logR_Copy_Number := logR_Copy_Number.ichor]
  titan_ichor_merged[is.na(Copy_Number), Copy_Number := Corrected_Copy_Number.ichor]
  titan_ichor_merged[is.na(MinorCN), MinorCN := 0]
  titan_ichor_merged[is.na(MajorCN), MajorCN := Corrected_Copy_Number.ichor]
  
  ## Fixing Other Parameters
  titan_ichor_merged[is.na(Median_logR), Median_logR := Median_logR.ichor] # assign ichor Median_logR to missing titan SNPs
  # titan_ichor_merged[is.na(Length.snp.), Length.snp. := Length.snp..ichor] # assign ichor Length.snp to missing titan SNPs - Good Idea?
  titan_ichor_merged[is.na(TITAN_call), TITAN_call := TITAN_call.ichor]
  titan_ichor_merged[is.na(Corrected_logR), Corrected_logR := Corrected_logR.ichor]
  
  ## Remove ichor tags
  titan_ichor_merged <- titan_ichor_merged[, -c(grep("ichor", colnames(titan_ichor_merged),value=T)), with=F]
  
  titan_ichor_merged$TITAN_call[titan_ichor_merged$Corrected_Call == 'HETD'] <- 'DLOH'
  titan_ichor_merged$TITAN_state[titan_ichor_merged$TITAN_call == 'DLOH'] <- 1
  titan_ichor_merged$TITAN_state[titan_ichor_merged$TITAN_call == 'HET'] <- 3
  titan_ichor_merged$TITAN_state[titan_ichor_merged$TITAN_call == 'GAIN'] <- 5
  
  segs2 <- extendSegments(titan_ichor_merged, removeCentromeres = TRUE, centromeres = centromeres, extendToTelomeres = FALSE,
                          chrs = chrs, genomeStyle = genomeStyle)
  
  
  ##################################### FIXING CN FILE OUTPUT #####################################
  correctCN2 <- correctIntegerCN(cn, segs2, purity, ploidyT, maxCNtoCorrect.autosomes = NULL, 
                                 maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.05, gender = gender, chrs = chrs)
  cn2 <- correctCN2$cn
  cn3 <- fill_cn_blanks_if_ichor(cn2)
  segs3 <- correctCN2$segs
  
  ##################################### Format for PLOTTING #####################################
  cn_plotting <- subset(cn3, ((!is.na(cn3[,"Position"])) & (!is.na(cn3[,"logR_Copy_Number"]))))
  cn_plotting <- cn_plotting[,c("Chr", "Position", "RefCount", "Depth", "AllelicRatio", "LogRatio", "CopyNumber", "TITANstate", "TITANcall", "ClonalCluster", "CellularPrevalence", "logR_Copy_Number", "Corrected_logR", "Corrected_Ratio", "Corrected_Copy_Number", "Corrected_Call")]
  
  segs_plotting <- copy(segs3)
  setnames(segs_plotting, c("Start", "End"), c("Start_Position.bp.", "End_Position.bp."))
  
  cn_plotting$Chr <- factor(cn_plotting$Chr, levels = chrs)
  cn_plotting <- cn_plotting[order(cn_plotting$Chr), ]
  
  ## Create a new column for segment-based copy number in the cn_plotting category
  seg_copy_numbers <- numeric(length = nrow(cn_plotting))
  for (i in 1:nrow(cn_plotting)) {
    # Check if POSITION falls between START and END in table 2
    matched_row <- which(cn_plotting$Chr[i] == segs_plotting$Chromosome & cn_plotting$Position[i] >= segs_plotting$Start_Position.bp. & cn_plotting$Position[i] <= segs_plotting$End_Position.bp.)
    # If a match is found, extract the copy number
    if (length(matched_row) > 0) {
      seg_copy_numbers[i] <- segs_plotting$Corrected_Copy_Number[matched_row]
    } else {
      seg_copy_numbers[i] <- cn_plotting$Corrected_Copy_Number[i]  # If no match found, store NA
    }
  }
  cn_plotting$Segment_Copy_Number <- seg_copy_numbers
  
  ## Load titanCNA Rdata
  load(titan_rdata)
  
  ## Plotting Genome-Wide Titan-Ichor Repaired File 
  numClustersStr <- as.character(numClusters)
  norm <- tail(convergeParams$n,1)
  chrs_plot <- c(chrs, "chrX")


  ## Plot Segment Function
  plotCNlogRByChr_segment <- function(dataIn, chr = c(1:22), segs = NULL, 
                                      plotCorrectedCN = TRUE, geneAnnot = NULL,
                                      ploidy = NULL, normal = NULL, spacing = 4, alphaVal = 1, xlim = NULL, ...) {
    # color coding
    alphaVal <- ceiling(alphaVal * 255)
    class(alphaVal) = "hexmode"
    cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",
                        "#BB0000", "#CC0000", "#DD0000", "#EE0000", rep("#FF0000",493))
                        cnCol <- paste(cnCol, alphaVal, sep = "")
                        # cnCol <-
                        # col2rgb(c('green','darkgreen','blue','darkred','red','brightred'))
                        names(cnCol) <- c(0:500)
                        
                        # use consistent chromosome naming convention
                        chr <- as.character(chr)
                        genomeStyle <- seqlevelsStyle(as.character(dataIn$Chr))[1]
                        chr <- mapSeqlevels(chr, genomeStyle, drop = FALSE)[1, ]
                        
                        if (plotCorrectedCN && "Corrected_Copy_Number" %in% colnames(dataIn)){
                          binCN <- "Segment_Copy_Number"
                          segCN <- "Segment_Copy_Number"
                        }else{
                          binCN <- "CopyNumber"
                          segCN <- "Copy_Number"
                        }
                        
                        dataIn <- data.table(dataIn)
                        ## adjust logR values for ploidy ##
                        if (!is.null(ploidy)) {
                          if (is.null(normal)){
                            stop("plotCNlogRByChr: Please provide \"normal\" contamination estimate.")
                          }
                          dataIn[, LogRatio := LogRatio + log2(((1-normal)*ploidy+normal*2)/2)]
                          
                          if (!is.null(segs)){
                            segs.sample <- copy(segs)
                            segs.sample[, Median_logR := Median_logR + log2(((1-normal)*ploidy+normal*2) / 2)]
                          }
                      }
                        # plot for all chromosomes specified
                        dataIn <- dataIn[Chr %in% chr, ]
                        coord <- getGenomeWidePositions(dataIn[, Chr], dataIn[, Position]) # Original line
                        #coordEnd <- getGenomeWidePositions(dataIn[,Chr],dataIn[, End]) # Added line
                        plot(coord$posns, as.numeric(dataIn[, logR_Copy_Number]),
                             col = cnCol[as.character(dataIn[, get(binCN)])],
                             pch = 16, xaxt = "n", las = 1, bty = "n",
                             ylab = "Copy Number (Segment)", ...)
                        lines(as.numeric(c(1, coord$posns[length(coord$posns)])),
                              rep(0, 2), type = "l", col = "grey", lwd = 2)
                        plotChrLines(dataIn[, Chr], coord$chrBkpt, par("yaxp")[1:2])
                        #plot segments
                        if (!is.null(segs)){
                          segs.sample <- segs.sample[Chromosome %in% chr]
                          coordEnd <- getGenomeWidePositions(segs.sample[, Chromosome], segs.sample[, End_Position.bp.])
                          coordStart <- coordEnd$posns - (segs.sample[, End_Position.bp.] - segs.sample[, Start_Position.bp.] + 1)
                          xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
                          col <- cnCol[as.character(segs.sample[, get(segCN)])]
                          value <- as.numeric(segs.sample[, Median_logR])
                          mat <- as.data.frame(cbind(coordStart, coordEnd$posns, value, col))
                          rownames(mat) <- 1:nrow(mat)
                          tmp <- apply(mat, 1, function(x){
                            lines(x[1:2], rep(x[3], 2), col = x[4], lwd = 3, lend = 1)
                          })
                        }
                      }

  # New limits for the logR Copy Number plotting
  plotYlim_CN <- c(0, 12)
  outFile <- paste0(outplot_repaired, id, "_cluster", numClustersStr, "_TitanIchor_Repaired_logRCN.pdf")
  #png(outFile,width=1000,height=300)
  pdf(paste0(outFile),width=20,height=6)
  plotCNlogRByChr_segment(dataIn=cn_plotting, chr=chrs_plot, segs = NULL, ploidy=ploidy,
                          normal = norm, geneAnnot=NULL, spacing=4, main=id, xlab="",
                          ylim=plotYlim_CN, cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, plotCorrectedCN = TRUE)
  dev.off()
  
  segs_out <- copy(segs3)
  cn_out <- copy(cn3)
  
  write.table(cn_plotting, file = paste0(outplot_repaired, id, "_cluster", numClustersStr, ".TitanIchor.Repaired.cna.plotting.txt"), col.names=T, row.names=F, quote=F, sep="\t")
  write.table(segs_plotting, file = paste0(outplot_repaired, id, "_cluster", numClustersStr, ".TitanIchor.Repaired.segs.plotting.txt"), col.names=T, row.names=F, quote=F, sep="\t")
  } else {
  ##################################### Runs if opted for no ichor rescue #####################################
  ## extend segments to remove gaps
  segs_out <- extendSegments(segs, removeCentromeres = TRUE, centromeres = centromeres, extendToTelomeres = FALSE,
                          chrs = chrs, genomeStyle = genomeStyle)
  cn_out <- copy(cn)
}

## write segments to file ##
write.table(segs_out, file = outSegFile, col.names=T, row.names=F, quote=F, sep="\t")
write.table(cn_out, file = outBinFile, col.names=T, row.names=F, quote=F, sep="\t")
## write segments without germline SNPs
outSegNoSNPFile <- gsub(".txt", ".noSNPs.txt", outSegFile)
write.table(segs3[, -c("Start.snp", "End.snp")], file = outSegNoSNPFile, col.names=T, row.names=F, quote=F, sep="\t")

## write segments in IGV / GISTIC format ##
chrs <- c(1:22, "X")
igv <- segs3[, .(Sample, Chromosome, Start.snp, End.snp, Length.snp., logR_Copy_Number)]
igv[Chromosome %in% chrs[1:22], Corrected.logR := log2(logR_Copy_Number / 2)]
igv[Chromosome == chrs[grep("X", chrs)], Corrected.logR := log2(logR_Copy_Number / 1)]
igv[, logR_Copy_Number := NULL]
outIGVFile <- gsub("seg.txt", "segIGV.txt", outSegFile)
write.table(igv, file = outIGVFile, col.names=T, row.names=F, quote=F, sep="\t")

save.image(outImageFile)


