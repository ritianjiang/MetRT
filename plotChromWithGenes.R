#' Plot chromosome and genes
#'
#' @export

plotChromosomeWithGenes <- function(chrom,ideo,bpLim=NULL,
                           titleVertical=1,titleHorizontal=4,vertical=TRUE,addScale=TRUE,
                           chromName_cex=0.8,verbose=FALSE,GeneList = NULL,...)
{
  ideo<-ideo[ideo$chrom == chrom,]
  chromStart = min(ideo$chromStart)
  chromEnd = max(ideo$chromEnd)
  chrom_width <- c(-9, 9)
  bpLim = c(chromStart,chromEnd)
  if (verbose) cat("\n plotting...",  chrom, "\n")
  args <- list(type='n',
               yaxt='n',xaxt='n',xlab=' ', ylab=' ',
               bty='n',x=0,y=0,
               cex.axis=1,font.axis=2,cex.lab=1.3,font=2)

  if (vertical) {
    args <- modifyList(args, list(yaxs='i', ylim=rev(bpLim),
                                  xlim=chrom_width))
    #mar=replace(def_mar,1,0.1)))
  } else {
    args <- modifyList(args, list(xaxs='i', xlim=bpLim,
                                  ylim=chrom_width))
    #mar=replace(def_mar,1,0.1)))
  }
  #browser()
  suppressWarnings(do.call('plot',args))
  if (addScale) {
    if(vertical) {
      at <-2
    } else {
      at <- 1
    }
    axt <- axTicks(side=at)
    suppressWarnings(axis(at, xpd=FALSE,at=axt, labels=axt/1000000, ...))
  }

  #drawing centromere
  centro_idx <- which(ideo$gieStain=="acen") # gieStain values of "acen"
  # represent the centromere
  if(any(centro_idx)){
    centro_coord <- range(c(ideo[centro_idx, "chromStart"],
                            ideo[centro_idx, "chromEnd"]))
    ideo <- ideo[-centro_idx, ]

    if(vertical) polygon(c(0, 1, 0, 1), c(centro_coord[1],
                                          centro_coord[2], centro_coord[2], centro_coord[1]),
                         col="darkred", border="black", lwd=2)
    else polygon(c(centro_coord[1], centro_coord[2], centro_coord[2],
                   centro_coord[1]), c(0, 1, 0, 1), col="darkred",
                 border="black", lwd=2)
    rm(centro_coord)

  };
  rm(centro_idx)

  ##setting band colors
  ## gneg & gposN(where N is an integer [1,100] ) prepresent densities.
  ## gneg has density 0
  col <- character(nrow(ideo))
  col[1:length(col)] <- "#000000" # gvar(constitutive heterochromatins)
  # and gpos100
  # colors will become darker as density increases from 0 -> 100
  col[grep("gneg", ideo$gieStain)] <- "gray90"
  col[grep("gpos25", ideo$gieStain)] <- "gray75"
  col[grep("gpos33", ideo$gieStain)] <- "gray66"
  col[grep("gpos50", ideo$gieStain)] <- "gray50"
  col[grep("stalk", ideo$gieStain)] <- "darkred" # repetitive areas
  col[grep("gpos66", ideo$gieStain)] <- "gray33"
  col[grep("gpos75", ideo$gieStain)] <- "gray25"
  ##--------------------------------------------------
  #adding bands
  for (k in 1:nrow(ideo)) {
    if(vertical) {
      rect(xleft=0, ybottom=ideo[k,2],xright=1, ytop=ideo[k,3],
           col=col[k], border=col[k])
    }
    else {
      rect(xleft=ideo[k,2],ybottom=0,xright=ideo[k,3],ytop=1,
           col=col[k], border=col[k])
    }

  }; rm(col)
  # add chrom name
  if (vertical) {
    mtext(chrom, side=titleVertical,line=2, outer=FALSE,
          col='gray50',adj=c(1,1),cex=chromName_cex,font=2)
  } else {
    mtext(chrom, side=titleHorizontal,line=0.5,outer=FALSE,
          col='gray50', cex=chromName_cex,font=2,las=3)
  }

  #add border
  p_arm <- grep("p", ideo$name)
  if(any(p_arm)){
    p_rect_coord <- range(c(ideo[p_arm, "chromStart"],
                            ideo[p_arm, "chromEnd"]))
    if(vertical) rect(0, p_rect_coord[1], 1, p_rect_coord[2],
                      col=NA, border="black", lwd=2)
    else rect(p_rect_coord[1], 0, p_rect_coord[2], 1,
              col=NA, border="black", lwd=2)
    rm(p_rect_coord)
  }
  q_arm <- grep("q", ideo$name)
  if(any(q_arm)){
    q_rect_coord <- range(c(ideo[q_arm, "chromStart"],
                            ideo[q_arm, "chromEnd"]))
    if(vertical) rect(0, q_rect_coord[1], 1, q_rect_coord[2],
                      col=NA, border="black", lwd=2)
    else rect(q_rect_coord[1], 0, q_rect_coord[2], 1, col=NA ,
              border="black", lwd=2)
    rm(q_rect_coord)
  }

  lableMin<-min(GeneList$start)-2000000
  labelMax<-max(GeneList$start)+2000000
  GeneList<-GeneList[order(GeneList$start),]
  step<-(labelMax - lableMin)/nrow(GeneList)
  for(g in 1:nrow(GeneList)){
    text(GeneList[g,]$gene,
         x = lableMin+(g-1)*step,y=6.5,srt=90,cex=0.8)
    segments(x0 = GeneList[g,]$start,y0 = 0,
             x1 = GeneList[g,]$start,y1 = 3)
    segments(x0 = GeneList[g,]$start,y0 = 3,
             x1 = lableMin + (g-1)*step,y1 = 5)
    # if(g%%2 == 1){
    #     if(g==1 | g%%4 ==1){
    #         segments(x0 = GeneList[g,]$start,y0 = 1,
    #                  x1 = GeneList[g,]$start + 1000,y1 = 5)
    #         text(GeneList[g,]$gene,
    #              x = GeneList[g,]$start+1000, y = 6.5,srt = 90,cex = 0.8)
    #     }
    #     else{segments(x0 = GeneList[g,]$start,y0 = 1,
    #                   x1 = GeneList[g,]$start + 1000,y1 = 6)
    #          text(GeneList[g,]$gene,
    #                   x = GeneList[g,]$start+1000, y = 7.5,srt = 90,cex = 0.8)}
    #   }
    # else{
    #   if(g==1 | g%%4 ==2){
    #     segments(x0 = GeneList[g,]$start,y0 = 0,
    #              x1 = GeneList[g,]$start + 1000,y1 = -4)
    #     text(GeneList[g,]$gene,
    #          x = GeneList[g,]$start+1000, y = -5.5,srt = 90,cex = 0.8)
    #   }
    #   else{segments(x0 = GeneList[g,]$start,y0 = 0,
    #                 x1 = GeneList[g,]$start + 1000,y1 = -6)
    #     text(GeneList[g,]$gene,
    #          x = GeneList[g,]$start+1000, y = -7.5,srt = 90,cex = 0.8)}
    # }


  }
}
