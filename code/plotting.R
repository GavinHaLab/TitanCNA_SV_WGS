###########################################
############ PLOTTING FUNCTION ###########
###########################################
## sv is a data.table output by combineSVABAandTITAN.R
### TO-DO: 
# 1) Filter by support
# 2) Strand direction of breakpoint
plotRearrangementArcs <- function (sv, cn, ploidy = NULL, interchr=TRUE, xlim=NULL, segIn=NULL, support=1, chr=NULL, chrLens=NULL, minSPAN = 10, buffer=1e6, plotAtCentre = FALSE, centreLine=0, arcHeight=4, lty=1, lcol = "black", svTypeCol=FALSE, arr.col = "black", arr.pos = 1, endhead = FALSE, lwd = 1, orient="topbottom", include.inter.chr = FALSE, offset.factor = 1.15){
	require(diagram)
	sv <- sv[sv$SPAN >= minSPAN | sv$SPAN == -1, ]
	sv <- as.data.frame(sv)
	
	# adjust for ploidy #
	if (!is.null(ploidy)){
    cn[, "LogRatio"] <- as.numeric(cn[, "LogRatio"]) + log2(ploidy / 2)
  }
	
  if (!is.null(chr) && nrow(sv) > 0){
    dataByChr <- sv[sv[,"chromosome_1"]==as.character(chr) | 
                  sv[,"chromosome_2"]==as.character(chr),]
    if (nrow(dataByChr) > 0){
      #use only interchromosomal rearrangements
      dataByChr$intraChr <- dataByChr[,"chromosome_1"] == dataByChr[,"chromosome_2"]
      dataByChr$inWindow <- "NONE"
      if (!is.null(xlim)){ 
        ## only include SV if at least one of the breakpoints is inside window
        leftIn <- (dataByChr[, "chromosome_1"] == chr & 
                      dataByChr[, "start_1"] >= xlim[1] - buffer & 
                      dataByChr[, "start_1"] <= xlim[2] + buffer)
        rightIn <- (dataByChr[, "chromosome_2"] == chr & 
                      dataByChr[, "start_2"] >= xlim[1] - buffer & 
                      dataByChr[, "start_2"] <= xlim[2] + buffer)
        dataByChr$inWindow[leftIn | rightIn] <- "ONE"
        ## SV for both breakpoints contained within window
        dataByChr$inWindow[leftIn & rightIn] <- "BOTH"						
      }		

		## plot intra-chromosomal SV ##
		
			for (i in 1:nrow(dataByChr)){
				if (dataByChr[i, "inWindow"] == "NONE"){
					next
				}
				#print(dataByChr[i, ])
				## get cn logR of nearest bin
				yLogR <- findNearestLogR(dataByChr[i, ], cn, buffer = buffer * 100)	
				yToUse <- yLogR			
				orient <- "top"
				start1 <- dataByChr[i,"start_1"]
        end1 <- dataByChr[i,"start_2"]

				if (plotAtCentre && !is.null(centreLine)){
					yToUse <- c(centreLine, centreLine)
				}
				curveHeight <- arcHeight[2] - mean(yToUse)#/ abs(start1 - end1)
				text.height <- arcHeight[2]
				if (is.null(arr.col)){
					arr.col.toUse <- arr.col.sv
				}else {arr.col.toUse <- arr.col}
				if (is.null(lcol)){
					lcol.toUse <- arr.col.sv
				}else {lcol.toUse <- lcol}
				if (svTypeCol){
					lcol.toUse <- dataByChr[i, "color"]
				}
				## if only 1 breakpoint in window ##
				if (dataByChr[i, "inWindow"] == "ONE"){
					## if inter-chromosomal ##
					if (dataByChr[i, "chromosome_1"] != chr || dataByChr[i, "chromosome_2"] != chr){
						orient <- "bottom"							
						text.height <- 1/arcHeight[2]			
						if (!is.null(centreLine)) { yToUse <- c(centreLine, centreLine) }
						curveHeight <- arcHeight[1] + mean(yToUse)
					}
	
				 	if (is.null(offset.factor)) { offset.factor <- 1.15 }
					leftBreak <- paste0(dataByChr[i,"chromosome_1"], ":", dataByChr[i, "start_1"])
					rightBreak <- paste0(dataByChr[i,"chromosome_2"], ":", dataByChr[i, "start_2"])
				 	## for inter chromosome
				 	interchr.offset <- (par("usr")[2] - par("usr")[1]) * offset.factor # left coordinate if outside plot
					if (dataByChr[i, "chromosome_1"] != chr){					    
						if (abs(dataByChr[i, "start_2"]-xlim[1]) < abs(dataByChr[i, "start_2"]-xlim[2])){
						# 2nd bkpt is closer to chr begin
							start1 <- par("usr")[1] - interchr.offset # left boundary  
							#text(x=par("usr")[1], y=text.height, pos=4, cex=0.5, label=leftBreak) 
						}else{
							start1 <- end1
							end1 <- par("usr")[2] + interchr.offset#right boundary  
							#text(x=par("usr")[2], y=text.height, pos=2, cex=0.5, label=rightBreak)
							#yToUse <- yToUse[c(2,1)]              
						}                   
					}else if (dataByChr[i, "chromosome_2"] != chr){
						if (abs(dataByChr[i, "start_1"]-xlim[1]) < abs(dataByChr[i, "start_1"]-xlim[2])){
						# 1st bkpt is closer to chr begin
							end1 <- start1
							start1 <- par("usr")[1] - interchr.offset #- (par("usr")[2] * offset.factor - par("usr")[2]) # left boundary 
							#text(x=par("usr")[1], y=text.height, pos=4, cex=0.5, label=leftBreak) 
							#yToUse <- yToUse[c(2,1)]               
						}else{
							end1 <- par("usr")[2] + interchr.offset #right boundary 
							#text(x=par("usr")[2], y=text.height, pos=2, cex=0.5, label=rightBreak)     
						}  
					}
					## for zoom within same chromosome
					if (dataByChr[i, "intraChr"]){
						if (start1 < xlim[1]){ 
							start1 <- par("usr")[1] - interchr.offset #- (par("usr")[2] * offset.factor - par("usr")[2]) # left boundary 
							#text(x=par("usr")[1], y=text.height, pos=4, cex=0.5, label=leftBreak) 
						}
						if (start1 > xlim[2]){
							start1 <- par("usr")[2] - interchr.offset #- (par("usr")[2] * offset.factor - par("usr")[2]) #right boundary
						}
						if (end1 < xlim[1]){
							end1 <- par("usr")[1] + interchr.offset # left boundary
						}
						if (end1 > xlim[2]){
							end1 <- par("usr")[2] + interchr.offset # left boundary   
							#text(x=par("usr")[2], y=text.height, pos=2, cex=0.5, label=rightBreak)             
						}        
						#if (sv.type == "inv") { end2 <- start1; start2 <- end1 }    
					}
					#curveHeight <- par("usr")[3]# - mean(yToUse)#* 2 / (end1 - start1)
				}
				
				curvedarrow.new(from = c(start1, yToUse[1]), endhead = endhead, orient = orient,
									to = c(end1, yToUse[2]), lty = lty, lwd = lwd, lcol=lcol.toUse,
									arr.col=arr.col.toUse, arr.width = 0.1, arr.length=0.2,
									curve = curveHeight, arr.type = "triangle", arr.pos = arr.pos)

			}
		}		
	}
	invisible()
}

## format of dataIn is output from /home/unix/gavinha/software/code/git/scripts/titan/analysis/combineTITAN-ichor.R
plotTitanIchorCNA <- function(dataIn, param = NULL, colName = "LogRatio", callColName="Corrected_Call", segs=NULL, chr=NULL, purity = NULL, ploidyT = NULL, geneAnnot=NULL, yrange=c(-4,6), yaxis = "logRatio", xlim=NULL, xaxt = "n", cex = 0.5, gene.cex = 0.5, plot.title = NULL, cnCol = NULL, spacing=4, cytoBand=T, genomewide=F, alphaVal=1, main){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  subcloneCol <- c("#00FF00")
  if (is.null(cnCol) & (genomewide==FALSE)){
    cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
    cnCol <- c(cnCol, "HET"="#0000FF", "DLOH"="#006400", "NLOH"="#0000FF", "ALOH"="#FF0000", "ASCNA"="#FF0000", "BCNA"="#FF0000", "UBCNA"="#FF0000")
    names(cnCol)[1:30] <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  } else if (is.null(cnCol) & (genomewide==TRUE) & ((callColName == "Segment_Copy_Number") | callColName == "Corrected_Copy_Number" | callColName == "Adjusted_Segment_Copy_Number")) {
    cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",rep("#FF0000", 26))
    cnCol <- paste(cnCol, alphaVal, sep = "")
    names(cnCol)[1:30] <- c(0,1,2,3,4,5,6:25)
  } else if (is.null(cnCol) & (genomewide == TRUE)) {
    cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
    cnCol <- c(cnCol, "HET"="#0000FF", "DLOH"="#006400", "NLOH"="#0000FF", "ALOH"="#FF0000", "ASCNA"="#FF0000", "BCNA"="#FF0000", "UBCNA"="#FF0000")
    names(cnCol)[1:30] <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
	} else {
		cnCol.col <- as.character(cnCol[1])
		cnCol <- c(cnCol, "HET"=cnCol.col, "DLOH"=cnCol.col, "NLOH"=cnCol.col, "ALOH"=cnCol.col, "ASCNA"=cnCol.col, "BCNA"=cnCol.col, "UBCNA"=cnCol.col)
		names(cnCol)[1:30] <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
	}
  # adjust for ploidy #
  normCN <- 2
  if (!is.null(ploidyT) & yaxis != "integer"){
    ploidyS <- purity * ploidyT + (1-purity) * normCN
    dataIn[, colName] <- as.numeric(dataIn[, colName]) + log2(ploidyS / 2)
    
    if (!is.null(segs)){
      segs[, colName] <- segs[, colName] + log2(ploidyS / 2)
    }
  }
  
  if (!is.null(chr)){
    if (genomewide) {
      dataByChr <- dataIn[dataIn[,"Chr"] %in% chr, ]
      dataByChr <- subset(dataByChr, !is.na(dataByChr[, colName]))
       
      zero <- 0.5  
      cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      #ploidyToUse <- ploidyS   
      if (yaxis == "integer"){
        y.ticks <- log2(cn)
        y.ticks[1] <- log2(zero)  
        yrange[1] <- y.ticks[1]    
        ylab <- "Copy Number (logR)"
        #dataByChr[, colName] <- log2(logRbasedCN(dataByChr[, colName], purity, ploidyT, cn=normCN))
        dataByChr[, colName] <- log2(dataByChr[, colName])
        if (colName != "LogRatio") {    ## Added conditional for this line
          dataByChr[dataByChr[[colName]] < -1 , colName] <- -1
        }
        # dataByChr[dataByChr[[colName]] < -1 , colName] <- -1
        if (!is.null(segs)){
      		segs[, colName] <- log2(segs[, colName])# + log2(ploidyS / 2)
    		}
        centreLine <- log2(normCN)
      }else{      
      	#dataByChr[, colName] <- dataByChr[, colName] + log2(ploidyS / 2)
      	cnLog <- log2(cn[-which(cn==3)] / normCN)  
        cn <- seq(-2,yrange[2],2)#c(-2, cn)
        y.ticks <- cn
        ylab <- "Copy Number (log2 ratio)"
        centreLine <- 0
      }
      # plot the data
      par(mar=c(spacing,8,4,2))
      dataByChr_table <- as.data.table(dataByChr)
      coord <- getGenomeWidePositions(dataByChr_table[, Chr], dataByChr_table[, Position]) # Original line

      if (callColName == "Segment_Copy_Number") {
        dataByChr_table[, Segment_Copy_Number := Segment_Copy_Number + 1]
        # dataByChr_table[grepl("X",Chr), Segment_Copy_Number := Segment_Copy_Number + 2]
        dataByChr <- as.data.frame(dataByChr_table)
      } else if (callColName == "Corrected_Copy_Number") {
          dataByChr_table[, Corrected_Copy_Number := Corrected_Copy_Number + 1]
          dataByChr <- as.data.frame(dataByChr_table)
      } else if (callColName == "Adjusted_Segment_Copy_Number") {
          dataByChr_table[, Adjusted_Segment_Copy_Number := Adjusted_Segment_Copy_Number + 1]
          dataByChr <- as.data.frame(dataByChr_table)
      } else {
        dataByChr <- as.data.frame(dataByChr_table)
      }
      col <- cnCol[dataByChr[,callColName]]
      plot(coord$posns, as.numeric(dataByChr[, colName]),
        col=col,
        pch=16, ylim=yrange, yaxt="n",
        xlim=xlim, xaxt = xaxt, xlab="",ylab=ylab,
        cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
      axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)

      lines(as.numeric(c(1, coord$posns[length(coord$posns)])), rep(centreLine,2), type = "l", col = "grey", lwd = 2)

      plotChrLines(dataByChr_table[, Chr], coord$chrBkpt, par("yaxp")[1:2])
    
    } else {
      for (i in chr){
        dataByChr <- dataIn[dataIn[,"Chr"]==as.character(i),]
        ## set y axis labels as either integer or logR copy number
        #avgTumPloidy <- round(ploidyT)
  
        zero <- 0.5  
        cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
        #ploidyToUse <- ploidyS
        if (i == "X"){
          normCN <- 1
          zero <- 0.25
          cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
        }      
        if (yaxis == "integer"){
          y.ticks <- log2(cn)
          y.ticks[1] <- log2(zero)  
          yrange[1] <- y.ticks[1]    
          ylab <- "Copy Number"
          #dataByChr[, colName] <- log2(logRbasedCN(dataByChr[, colName], purity, ploidyT, cn=normCN))
          dataByChr[, colName] <- log2(dataByChr[, colName])
          dataByChr[dataByChr[[colName]] < -1 , colName] <- -1
          if (!is.null(segs)){
            segs[, colName] <- log2(segs[, colName])# + log2(ploidyS / 2)
          }
          centreLine <- log2(normCN)
        }else{      
          #dataByChr[, colName] <- dataByChr[, colName] + log2(ploidyS / 2)
          cnLog <- log2(cn[-which(cn==3)] / normCN)  
          cn <- seq(-2,yrange[2],2)#c(-2, cn)
          y.ticks <- cn
          ylab <- "Copy Number (log2 ratio)"
          centreLine <- 0
        }

        #plot the data
        #if (outfile!=""){ pdf(outfile,width=10,height=6) }
        par(mar=c(spacing,8,4,2))
        #par(xpd=NA)
        coord <- (as.numeric(dataByChr[,"End"]) + as.numeric(dataByChr[,"Start"]))/2
        if (is.null(xlim)){
          xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Start"]))
          xaxt <- "n"
        }
        if (is.null(plot.title)){
          plot.title <- paste("Chromosome ",i,sep="")
        }
        ## plot logR for bins ##
        plot(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,callColName]],
            pch=16, ylim=yrange, yaxt="n",
            xlim=xlim, xaxt = xaxt, xlab="",ylab=ylab,
            cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
        axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
        title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
        ## plot centre line ##
        lines(c(1,tail(na.omit(dataByChr[,3]), 1)),rep(centreLine,2),type="l",col="grey",lwd=0.75)
        if (!is.null(segs)){
          segsByChr <- segs[segs[,"Chromosome"]==as.character(i),,drop=FALSE]
          #ind <- segsByChr$subclone.status == FALSE
          apply(segsByChr, 1, function(x){
            lines(x[c("Start","End")], rep(x[colName], 2), col = cnCol[x[callColName]], lwd = 3)
            invisible()
          })
          #if (sum(!ind) > 0){
          #  apply(segsByChr[!ind, ], 1, function(x){
          #    lines(x[c("Start","End")], rep(x["Median_logR"], 2), col = subcloneCol, lwd = 3)
          #    invisible()
          #  })
          #}
        }
      }
      
      if (cytoBand==TRUE){
        require(quantsmooth)
        par(xpd = NA)
        #paintCytobands(chrom=chr, units="bases", pos=c(0,(yrange[1]-0.5)), width=0.75, legend=F)	
      }
      
      if (!is.null(geneAnnot)){
        #par(xpd=F)
        colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
        geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
        if (nrow(geneAnnot) > 0){
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)			
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"Start"]){ atP <- dataByChr[1,"Start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"Start"]){ atP <- dataByChr[dim(dataByChr)[1],"Start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)
          }
        }
      }
    }
  }
}

## compute copy number using corrected log ratio ##
logRbasedCN <- function(x, purity, ploidyT, cn = 2){
	ct <- (2^x * (cn * (1 - purity) + purity * ploidyT * (cn / 2)) - cn * (1 - purity)) / purity
	ct <- sapply(ct, max, 1/2^6)
	return(ct)
}


## modified to work for combine_TITAN_ICHOR/titan_ichor_cn.txt
findNearestLogR <- function(x, y, buffer = 1e6){
	y1 <- 0; y2 <- 0
	#y <- na.omit(y)
	bin1Ind <- x[["chromosome_1"]] == y[, "Chr"] & 
							x[["start_1"]] >= (y[, "Start"]) & 
							x[["start_1"]] <= (y[, "End"] + buffer)
	bin2Ind <- x[["chromosome_2"]] == y[, "Chr"] & 
							x[["start_2"]] >= (y[, "Start"]) & 
							x[["start_2"]] <= (y[, "End"] + buffer)
	if (sum(bin1Ind, na.rm=T) > 0){
	  # look at closest left and right points, find minimum distance to points, assign copy of min point
	  #ind <- (tail(which(bin1Ind), 1)-1):(tail(which(bin1Ind), 1))
	  ind <- tail(which(bin1Ind), 1)
	  #indClose <- ind[which.min(c(abs(y[ind[1], "end"]-x[["start_1"]]), abs(y[ind[2], "start"]-x[["start_1"]])))]
	  #y1 <- y[indClose, "copy"]
		#y1 <- max(y[tail(which(bin1Ind), 1):(tail(which(bin1Ind), 1)+1), "copy"], na.rm = TRUE)
	if (x[["orient_1"]] == "rev"){ # region to left of breakpoint 
      y1 <- y[max(ind - 1, 1), "LogRatio"]
      #y1 <- tail(na.omit(y[max((ind-10):ind), "LogRatio"]), 1)
    }else if (x[["orient_1"]] == "fwd"){  # region to right of breakpoint 
      y1 <- y[min(ind + 1, length(bin1Ind)), "LogRatio"]
      #y1 <- head(na.omit(y[max(ind:(ind+10)), "LogRatio"]), 1)
    }
	}
	if (sum(bin2Ind, na.rm=T) > 0){
	  #ind <- (tail(which(bin2Ind), 1)-1):(tail(which(bin2Ind), 1))
	  ind <- tail(which(bin2Ind), 1)
	  #indClose <- ind[which.min(c(abs(y[ind[1], "end"]-x[["start_1"]]), abs(y[ind[2], "start"]-x[["start_1"]])))]
	  #y2 <- y[indClose, "copy"]
		#y2 <- max(y[tail(which(bin2Ind), 1):(tail(which(bin2Ind), 1)+1), "copy"], na.rm = TRUE)
		if (x[["orient_2"]] == "rev"){  # region to left of breakpoint 
      y2 <- y[max(ind - 1, 1), "LogRatio"]
    }else if (x[["orient_2"]] == "fwd"){  # region to right of breakpoint 
      y2 <- y[min(ind + 1, length(bin1Ind)), "LogRatio"]
    }
	}
	if (is.na(y1)) { y1 <- 0 }
	if (is.na(y2)) { y2 <- 0 }
	return(c(y1, y2))
}

curvedarrow.new <- function (from, to, lwd = 2, lty = 1, lcol = "black", arr.col = lcol,
    arr.pos = 0.5, curve = 1, dr = 0.01, endhead = FALSE, segment = c(0, 1), orient = "top", ...)
{
    dpos <- to - from
    angle <- atan(dpos[2]/dpos[1]) * 180/pi
    if (is.nan(angle))
        return
    mid <- 0.5 * (to + from)
    dst <- dist(rbind(to, from))
    ry <- abs(curve) #* dst
    aFrom <- 0
    aTo <- pi
    if (orient == "top") {
      aFrom <- 2 * pi
      aTo <- pi
    }else if (orient == "bottom"){
    	aFrom <- pi
    	aTo <- 2 * pi
    }
    if (segment[1] != 0)
        From <- segment[1] * aTo + (1 - segment[1]) * aFrom
    else From <- aFrom
    if (segment[2] != 1)
        To <- segment[2] * aTo + (1 - segment[2]) * aFrom
    else To <- aTo
    if ((from[1] < to[1] && orient == "top") || (from[1] > to[1] && orient == "bottom")){ ## draw arrow from down to upstream
      meanpi <- arr.pos * aFrom + (1 - arr.pos) * aTo
    }else{# if (from[1] > to[1]){ ## draw arrow from upstream to down 
      meanpi <- arr.pos * aTo + (1 - arr.pos) * aFrom
    }
    #if (endhead)
        #To <- meanpi 
    plotellipse(rx = dst/2, ry = ry, mid = mid, angle = angle,
        from = From, to = To, lwd = lwd, lty = lty, lcol = lcol)
    ell <- getellipse(rx = dst/2, ry = ry, mid = mid, angle = angle,
        from = 1.001 * meanpi, to = 0.999 * meanpi, dr = 0.002)
    if (endhead){
      if ((from[1] < to[1] && orient == "top") || (from[1] > to[1] && orient == "bottom")){ ## draw arrow from down to upstream
        code <- 2  # arrow draw at x1, y1
      }else{
        code <- 1 #arrow drawn at x0, y0
      }
      Arrows(ell[1, 1], y0=ell[1, 2], ell[nrow(ell), 1], y1=ell[nrow(ell), 2], 
          code = code, lcol = lcol, arr.col = arr.col, ...)
    }
    curvedarrow <- c(ell[nrow(ell), 1], ell[nrow(ell), 2])
}

### modify SNPchip function "plotIdiogram"
plotIdiogram.hg38 <- function (chromosome, cytoband, seqinfo, cytoband.ycoords, xlim,
                               ylim = c(0, 2), new = TRUE, label.cytoband = TRUE, label.y = NULL,
                               srt, cex.axis = 1, outer = FALSE, taper = 0.15, verbose = FALSE,
                               unit = c("bp", "Mb"), is.lattice = FALSE, ...)
{
  def.par <- par(no.readonly = TRUE, mar = c(4.1, 0.1, 3.1,
                                             2.1))
  on.exit(def.par)
  if (is.lattice) {
    segments <- lsegments
    polygon <- lpolygon
  }
  
  cytoband <- cytoband[cytoband[, "chrom"] == chromosome, ]
  unit <- match.arg(unit)
  if (unit == "Mb") {
    cytoband$start <- cytoband$start/1e+06
    cytoband$end <- cytoband$end/1e+06
  }
  if (missing(cytoband.ycoords)) {
    cytoband.ycoords <- ylim
  }
  rownames(cytoband) <- as.character(cytoband[, "name"])
  sl <- seqlengths(seqinfo)[chromosome]
  if (missing(xlim))
    xlim <- c(0, sl)
  if (unit == "Mb")
    xlim <- xlim/1e+06
  cytoband_p <- cytoband[grep("^p", rownames(cytoband), value = TRUE),
                         ]
  cytoband_q <- cytoband[grep("^q", rownames(cytoband), value = TRUE),
                         ]
  p.bands <- nrow(cytoband_p)
  cut.left <- c()
  cut.right <- c()
  for (i in seq_len(nrow(cytoband))) {
    if (i == 1) {
      cut.left[i] <- TRUE
      cut.right[i] <- FALSE
    }
    else if (i == p.bands) {
      cut.left[i] <- FALSE
      cut.right[i] <- TRUE
    }
    else if (i == (p.bands + 1)) {
      cut.left[i] <- TRUE
      cut.right[i] <- FALSE
    }
    else if (i == nrow(cytoband)) {
      cut.left[i] <- FALSE
      cut.right[i] <- TRUE
    }
    else {
      cut.left[i] <- FALSE
      cut.right[i] <- FALSE
    }
  }
  for (i in seq_len(nrow(cytoband))) {
    if (as.character(cytoband[i, "gieStain"]) == "stalk") {
      cut.right[i - 1] <- TRUE
      cut.left[i] <- NA
      cut.right[i] <- NA
      cut.left[i + 1] <- TRUE
    }
  }
  include <- cytoband[, "end"] > xlim[1] & cytoband[, "start"] <
    xlim[2]
  cytoband <- cytoband[include, ]
  N <- nrow(cytoband)
  cut.left <- cut.left[include]
  cut.right <- cut.right[include]
  if (new) {
    xx <- c(0, cytoband[nrow(cytoband), "end"])
    yy <- cytoband.ycoords
    plot(xx, yy, xlim = xlim, type = "n", xlab = "", ylab = "",
         axes = FALSE, yaxs = "i", ylim = ylim, ...)
  }
  top <- cytoband.ycoords[2]
  bot <- cytoband.ycoords[1]
  h <- top - bot
  p <- taper
  for (i in seq_len(nrow(cytoband))) {
    start <- cytoband[i, "start"]
    last <- cytoband[i, "end"]
    delta = (last - start)/4
    getStain <- function(stain) {
      switch(stain, gneg = "grey100", gpos25 = "grey90",
             gpos50 = "grey70", gpos75 = "grey40", gpos100 = "grey0",
             gvar = "grey100", stalk = "brown3", acen = "brown4",
             "white")
    }
    color <- getStain(as.character(cytoband[i, "gieStain"]))
    if (is.na(cut.left[i]) & is.na(cut.right[i])) {
      delta <- (last - start)/3
      segments(start + delta, cytoband.ycoords[1], start +
                 delta, cytoband.ycoords[2])
      segments(last - delta, cytoband.ycoords[1], last -
                 delta, cytoband.ycoords[2])
    }
    else if (cut.left[i] & cut.right[i]) {
      yy <- c(bot + p * h, bot, bot, bot + p * h, top -
                p * h, top, top, top - p * h)
      polygon(c(start, start + delta, last - delta, last,
                last, last - delta, start + delta, start), yy,
              col = color)
    }
    else if (cut.left[i]) {
      yy <- c(bot + p * h, bot, bot, top, top, top - p *
                h)
      polygon(c(start, start + delta, last, last, start +
                  delta, start), yy, col = color)
    }
    else if (cut.right[i]) {
      yy <- c(bot, bot, bot + p * h, top - p * h, top,
              top)
      polygon(c(start, last - delta, last, last, last -
                  delta, start), yy, col = color)
    }
    else {
      polygon(c(start, last, last, start), c(bot, bot,
                                             top, top), col = color)
    }
  }
  my.x <- (cytoband[, "start"] + cytoband[, "end"])/2
  if (label.cytoband & !is.lattice) {
    if (is.null(label.y)) {
      axis(1, at = my.x, labels = rownames(cytoband), outer = outer,
           cex.axis = cex.axis, line = 1, las = 3, tick = FALSE)
      axis(1, at = cytoband$start, outer = outer, cex.axis = cex.axis,
           line = 1, las = 3, labels = FALSE)
    }
    else {
      if (!is.numeric(label.y)) {
        warning("label.y must be numeric -- using default y coordinates for cytoband labels")
        label.y <- bot - p * h
      }
      if (missing(srt))
        srt <- 90
      text(x = my.x, y = rep(label.y, length(my.x)), labels = rownames(cytoband),
           srt = srt)
    }
  }
  return()
}


# Function to fill in CellularPrevalence for NA rows
fill_cellular_prevalence <- function(data, logRCN_col, CCN_col) {
  # Loop through dataframe
  for (i in 1:nrow(data)) {
    # Absolute restriction: Do not work with rows that have no value in ichor.logR_Copy_Number or are chrX
    if (data[i, "Chr"] == "chrX") {
      next
    } else {
      if ((is.na(data[i, ..logRCN_col])) & (i != 1)) {
        print(paste0("Skipping this row due to absent logR_Copy_Number: ", i))
        next
      }
      # Zeroth case: if the current row is the FIRST row, assign the CellularPrevalence value from the first non-NA row to this current row
      if (i == 1 & is.na(data[i, "CellularPrevalence"])) {
        first_non_na <- i
        while (is.na(data[first_non_na, "CellularPrevalence"])) {
          first_non_na <- first_non_na + 1
        }
        data[i, "CellularPrevalence"] <- data[first_non_na, "CellularPrevalence"]
      }
      # Same if the current row is the LAST row
      else if (i == nrow(data) & is.na(data[i, "CellularPrevalence"])) {
        last_non_na <- i
        while (is.na(data[last_non_na, "CellularPrevalence"])) {
          last_non_na <- last_non_na - 1
        }
        data[i, "CellularPrevalence"] <- data[last_non_na, "CellularPrevalence"]
      }
      # First Case: Where we can easily assign NA rows - if the current row has NA in CellularPrevalence and the value in the previous row is non-NA and is the SAME as the first non-NA row after current, assign the same CellularPrevalence to this row
      else if ((i != 1) & (i != nrow(data)) & (is.na(data[i, "CellularPrevalence"]))) {
        first_non_na_before <- i - 1
        while (is.na(data[first_non_na_before, "CellularPrevalence"])) {
          first_non_na_before <- first_non_na_before - 1
        }
        first_non_na_after <- i + 1
        while (is.na(data[first_non_na_after, "CellularPrevalence"])) {
          if ((first_non_na_after == nrow(data)) & (is.na(data[first_non_na_after, "CellularPrevalence"]))) {
            # Treat for the possibility that the last FEW rows are NA in CellularPrevalence; in that case, just fill in the rest with the same CellularPrevalence as the first_non_na_before
            data[first_non_na_after, "CellularPrevalence"] <- data[first_non_na_before, "CellularPrevalence"]
            break
          } 
          first_non_na_after <- first_non_na_after + 1
        }
        # if (first_non_na_after)
        if ((!is.na(data[i-1, "CellularPrevalence"]) & !is.na(data[i+1, "CellularPrevalence"]))) {
          if (data[i-1, "CellularPrevalence"] == data[first_non_na_after, "CellularPrevalence"]) {
            data[i, "CellularPrevalence"] <- data[i-1, "CellularPrevalence"]
            print(paste0("Row number: ", i, ", error is here in first case...............done"))
          # Second case: Where the current row has NA in CellularPrevalence and the value in the previous row is non-NA and the value is non-NA in the next row, and they are different. In this case, compare the id.logR_Copy_Number values, and assign the CellularPrevalence from the row with the id.logR_Copy_Number value closest to current row
          } else if ((data[i-1, "CellularPrevalence"] != data[i+1, "CellularPrevalence"])  & (!is.na(data[first_non_na_before, ..logRCN_col]) & (!is.na(data[first_non_na_after, ..logRCN_col])))) {
            if (!is.na(data[i, ..logRCN_col])) {
              diff_before <- abs(data[i, ..logRCN_col] - data[i-1, ..logRCN_col])
              diff_after <- abs(data[i, ..logRCN_col] - data[i+1, ..logRCN_col])
              if (diff_before < diff_after) {
                data[i, "CellularPrevalence"] <- data[i-1, "CellularPrevalence"]
              } else {
                data[i, "CellularPrevalence"] <- data[i+1, "CellularPrevalence"]
              }
            }
            print(paste0("Row number: ", i, ", error is here in second case...............done"))
          }
        }
        # Third-1 Case: if there are multiple NA rows in a row, we need to get the CellularPrevalence value for the first non-NA row before the NA rows, and if it matches the CellularPrevalence value to the first non-NA row after the NA rows, assign that same CellularPrevalence to all the NA rows
        else if ((!is.na(data[first_non_na_before, "CellularPrevalence"]) & is.na(data[i+1, "CellularPrevalence"]))) {
          if (data[first_non_na_before, "CellularPrevalence"] == data[first_non_na_after, "CellularPrevalence"]) {
            data[i, "CellularPrevalence"] <- data[first_non_na_before, "CellularPrevalence"]
            print(paste0("Row number: ", i, ", error is here in third-1 case...............done"))
          }
          # Third-2 case: if there are multiple NA rows, and the CellularPrevalence value is different for the first non-NA row before the NA rows compared to the CellularPrevalence value from the first non-NA row after the NA rows: take the average of the id.logR_Copy_Number values for the block of non-NA rows before current row that share the same CellularPrevalence up to the previous NA row or change in CellularPrevalence, and take the average of the id.logR_Copy_Number values for the block of non-NA rows after the current row that share the same CellularPrevalence, up to either the next NA row or change in CellularPrevalence. Compare the current row's id.logR_Copy_Number each of the two averages, and assign to this current row the CellularPrevalence value from the block of non-NA rows (before or after current row) that the current row's id.logR_Copy_Number is closest to.
          else if ((data[first_non_na_before, "CellularPrevalence"] != data[first_non_na_after, "CellularPrevalence"])) {
              # Get the average .logR_Copy_Number of all rows before current row that have the same CellularPrevalence, up until the row that either is NA in CellularPrevalence or has a different CellularPrevalence value
              check_cp_before <- TRUE
              if (!is.na(data[i, ..logRCN_col])) {
                sum_before <- data[first_non_na_before, ..logRCN_col]
                count_before <- 0
                latest_before_row <- first_non_na_before
                # Get the overall sum for .logR_Copy_Number and count of the rows that are above first_non_na_before that share the same CellularPrevalence value (and are non-NA)
                while (check_cp_before) {
                  if (latest_before_row != 1) {
                    before_row <- latest_before_row - 1
                    # print(paste0("........iterating row ", before_row))
                    if (!is.na(data[before_row, "CellularPrevalence"]) & (data[before_row, "CellularPrevalence"] == data[latest_before_row, "CellularPrevalence"])) {
                      sum_before <- sum_before + data[before_row, ..logRCN_col]
                      count_before <- count_before + 1
                      latest_before_row <- before_row
                    } else {
                      check_cp_before <- FALSE
                    }
                  }
                  check_cp_before <- FALSE
                }
                avg_before <- sum_before / count_before
                # Get the average id.logR_Copy_Number of all non-NA CellularPrevalence rows after current row, up until the next row that is NA in CellularPrevalence or has a different CellularPrevalence value
                check_cp_after <- TRUE
                sum_after <- data[first_non_na_after, ..logRCN_col]
                count_after <- 0
                latest_after_row <- first_non_na_after
                # Get the overall sum for id.logR_Copy_Number and count of the rows that are below after that share the same CellularPrevalence value (and are non-NA)
                while (check_cp_after) {
                  if (latest_after_row != nrow(data)) {
                    after_row <- latest_after_row + 1
                    if ((!is.na(data[after_row, "CellularPrevalence"])) & (data[after_row, "CellularPrevalence"] == data[latest_after_row, "CellularPrevalence"])) {
                      sum_after <- sum_after + data[after_row, ..logRCN_col]
                      count_after <- count_after + 1
                      latest_after_row <- after_row
                    } else {
                      check_cp_after <- FALSE
                    }
                  }
                  check_cp_after <- FALSE
                }
                avg_after <- sum_after / count_after
                diff_before <- abs(data[i, ..logRCN_col] - avg_before)
                diff_after <- abs(data[i, ..logRCN_col] - avg_after)
                if (diff_before < diff_after) {
                  data[i, "CellularPrevalence"] <- data[first_non_na_before, "CellularPrevalence"]
                } else {
                  data[i, "CellularPrevalence"] <- data[first_non_na_after, "CellularPrevalence"]
                }
              } else {
                  if (data[i, ..CCN_col] == data[first_non_na_before, ..CCN_col]) {
                    data[i, "CellularPrevalence"] <- data[first_non_na_before, "CellularPrevalence"]
                  } else if (data[i, ..CCN_col] == data[first_non_na_after, ..CCN_col]) {
                data[i, "CellularPrevalence"] <- data[first_non_na_after, "CellularPrevalence"]
                }
              }
              print(paste0("Row number: ", i, ", error is here in third-2 case...............done"))
          }
        } 
      } 
    } 
  }
  return (data)
}

# Function to adjust seg file positions based on ichorCNA insertion (from merging version)
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

# Function to fill in Position values for NA rows in bin file (that have data from ichorCNA)
fill_cn_blanks_if_ichor <- function(table) {
  table %>%
    mutate(Position = if_else(
      is.na(Position) & !is.na(.[[ncol(.)]]), 
      .$Start, 
      Position
    )) %>%
    return()
}


# Function to partition "empty" segments from the segment file, for processing with ichor bins
initialize_new_segments <- function(titan_segs_table, ichor_segs_table) {
  for (row in 1:(nrow(titan_segs_table) - 1)) {
  # For row 1 specifically: if the start of current row is greater than the End of row 1 in the cn2 file, then insert a row before the first row of all NA values except for the sample, chromosome, and start and end positions using start and end position from the cn2 file
    cur_row <- titan_segs_table[row,]
    first_ichor_seg_of_chr <- (ichor_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome,])[1,]
    last_ichor_seg_of_chr <- (ichor_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome,])[nrow(ichor_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome,]),]
    first_titan_seg_of_chr <- (titan_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome,])[1,]
    last_titan_seg_of_chr <- (titan_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome,])[nrow(titan_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome,]),]
    if (row == 1) {
      prev_row <- ichor_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome & Start < cur_row$Start,]
      if ((nrow(prev_row) != 0)) {
        # if (cur_row$Start > prev_row$End) {
        titan_segs_table <- rbind(data.table(Sample=cur_row$Sample, Chromosome=cur_row$Chromosome, Start=prev_row$Start, End=(cur_row$Start-1), Length.snp.=NA, Median_logR=NA, TITAN_call=NA, logR_Copy_Number=NA, Copy_Number=NA, Corrected_Copy_Number=NA, Corrected_Call=NA, Corrected_logR=NA, Corrected_MinorCN=NA, Corrected_MajorCN=NA, MinorCN=NA, MajorCN=NA, Sparse=1), titan_segs_table, fill=TRUE)
        # }
      }
  # If current row is the FIRST row in which the chromosome value changes to a new value, treat it the same case as row 1
    } else if ((cur_row$Start == first_titan_seg_of_chr$Start) & (cur_row$Start > first_ichor_seg_of_chr$Start) & (cur_row$End < first_ichor_seg_of_chr$End) & (row != 2)) {
    # if ((cur_row == first_titan_seg_of_chr) & (cur_row$Start > first_ichor_seg_of_chr$Start)) {
      prev_row <- ichor_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome & Start < cur_row$Start & End > cur_row$End,]
      if ((nrow(prev_row) != 0)) {
        # if (cur_row$Start > prev_row$End) {
          titan_segs_table <- rbind(titan_segs_table[1:(row-1),], data.table(Sample=cur_row$Sample, Chromosome=cur_row$Chromosome, Start=prev_row$Start, End=(cur_row$Start-1), Length.snp.=NA, Median_logR=NA, TITAN_call=NA, logR_Copy_Number=NA, Copy_Number=NA, Corrected_Copy_Number=NA, Corrected_Call=NA, Corrected_logR=NA, Corrected_MinorCN=NA, Corrected_MajorCN=NA, MinorCN=NA, MajorCN=NA, Sparse=1), titan_segs_table[row:nrow(titan_segs_table),], fill=TRUE)
        # }
      }
  # # If current row is the LAST row in which the chromosome value is the current row before itchanges to new row, treat similar to above case
  #   } else if ((cur_row$End == last_titan_seg_of_chr$End) & (cur_row$End < last_ichor_seg_of_chr$End)) {
  #     next_row <- ichor_segs_table[Sample == cur_row$Sample & Chromosome == cur_row$Chromosome & End > cur_row$End,]
  #     if ((nrow(next_row) != 0)) {
  #       # if (cur_row$End < next_row$Start) {
  #         titan_segs_table <- rbind(titan_segs_table[1:(row-1),], data.table(Sample=cur_row$Sample, Chromosome=cur_row$Chromosome, Start=(cur_row$End+1), End=next_row$End, Length.snp.=NA, Median_logR=NA, TITAN_call=NA, logR_Copy_Number=NA, Copy_Number=NA, Corrected_Copy_Number=NA, Corrected_Call=NA, Corrected_logR=NA, Corrected_MinorCN=NA, Corrected_MajorCN=NA, MinorCN=NA, MajorCN=NA, Sparse=1), titan_segs_table[row:nrow(titan_segs_table),], fill=TRUE)
  #       # }
  #     }
    } else {
      prev_row <- titan_segs_table[(row - 1),]
      next_row <- titan_segs_table[(row+1),]
      if ((all(sapply(list(prev_row$Chromosome, cur_row$Chromosome, next_row$Chromosome), function(x) x == cur_row$Chromosome))) & ((cur_row$Start - prev_row$End) > 100000)) {
        titan_segs_table <- rbind(titan_segs_table[1:(row-1),], data.table(Sample=cur_row$Sample, Chromosome=cur_row$Chromosome, Start=prev_row$End+1, End=cur_row$Start-1, Length.snp.=NA, Median_logR=NA, TITAN_call=NA, logR_Copy_Number=NA, Copy_Number=NA, Corrected_Copy_Number=NA, Corrected_Call=NA, Corrected_logR=NA, Corrected_MinorCN=NA, Corrected_MajorCN=NA, MinorCN=NA, MajorCN=NA, Sparse=1), titan_segs_table[row:nrow(titan_segs_table),], fill=TRUE)
      }
    }
  }
  return(titan_segs_table)
}


# Function to fill in the "empty" segments from the ichorCNA bin file
fill_new_segments <- function(titan_segs_table, bin_table) {
    for (row in 1:nrow(titan_segs_table)) {
    cur_row <- titan_segs_table[row,]
    if (cur_row$Sparse == 1) {
      bin_table_subset <- bin_table[which(cn2$Sample == cur_row$Sample & bin_table$Chr == cur_row$Chromosome & bin_table$Position >= cur_row$Start & bin_table$Position <= cur_row$End),]
      titan_segs_table[row, Length.snp. := sum(!is.na(bin_table_subset$Position))]
      # Create a further subset in which Position values are not NA
      bin_table_subset <- bin_table_subset[which(!is.na(bin_table_subset$Position)),]
      # Also take the median of their logR_Copy_Number values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, logR_Copy_Number := median(bin_table_subset$logR_Copy_Number, na.rm=TRUE)]
      # Also take the median of their Corrected_logR values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, Corrected_logR := median(bin_table_subset$Corrected_logR, na.rm=TRUE)]
      # Also take the median of their Corrected_Copy_Number values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, Corrected_Copy_Number := median(bin_table_subset$Corrected_Copy_Number, na.rm=TRUE)]
      # Also take the median of their Copy_Number values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, Copy_Number := median(bin_table_subset$CopyNumber, na.rm=TRUE)]
      # Also take the median of their Corrected_Call values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, Corrected_Call := median(bin_table_subset$Corrected_Call, na.rm=TRUE)]
      # Also take the median of their CellularPrevalence values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, Cellular_Prevalence := median(bin_table_subset$CellularPrevalence, na.rm=TRUE)]
      # Also take the median of their Corrected_Ratio values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, Corrected_Ratio := median(bin_table_subset$Corrected_Ratio, na.rm=TRUE)]
      # Also take th median of their LogRatio values for all rows in which Position is not NA and assign it to the sparse row
      titan_segs_table[row, Median_logR := median(bin_table_subset$LogRatio, na.rm=TRUE)]
    }
  }
  return(titan_segs_table)
}


# Function to check segment copy numbers in bin file, and correct if necessary
correct_segCN_bins <- function(bin_table, CorrectedCN_Col, SegCN_Col, threshold) {
  for (i in 1:nrow(bin_table)) {
    # Do the same with Segment_Copy_Number
    if (abs(bin_table$Segment_Copy_Number[i] - bin_table$Corrected_Copy_Number[i]) >= 2) {
      bin_table$Segment_Copy_Number[i] <- bin_table$Corrected_Copy_Number[i]
    }
    # Other condition: keep a counter of the number of rows in a row in which Segment_Copy_Number or NewMethod_Segment_Copy_Number differs from Corrected_Copy_Number; if that number is greater than 10, then replace Segment_Copy_Number or NewMethod_Segment_Copy_Number with Corrected_Copy_Number for all rows within the counter. Counter resets as soon as there is a row in which Segment_Copy_Number or NewMethod_Segment_Copy_Number is equal to Corrected_Copy_Number
    if (bin_table$Segment_Copy_Number[i] != bin_table$Corrected_Copy_Number[i]) {
      counter <- 1
      for (j in (i+1):nrow(cn_plotting)) {
        if (bin_table$Segment_Copy_Number[j] != bin_table$Corrected_Copy_Number[j]) {
          counter <- counter + 1
        } else {
          break
        }
      }
      if (counter > threshold) {
        for (j in i:(i+counter)) {
          bin_table$Segment_Copy_Number[j] <- bin_table$Corrected_Copy_Number[j]
        }
      }
    }
  }
  return(bin_table)
}