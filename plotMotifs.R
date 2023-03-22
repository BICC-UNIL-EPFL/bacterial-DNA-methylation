#setwd('E:/VMShared/MB/AB_A118/')
# Generate heatmap and circos plots from the detected nucleotide modification data
# Copyright (c) 2023 Christian Iseli, UNIL - EPFL Bioinformatics Competence Center
library(circlize)
library(gplots)
options("scipen"=6)
mycex <- 0.3
mycexlegend <- 0.5
mynotecex <- 0.5
axisFontSize <- 0.4
labelXpos <- 0.5
labelYpos <- 5.8
pad <- 2000
perpage <- 3
nbdiv <- perpage + 2

x <- read.table('all_motifs.txt',header=T)
a.all <- read.csv('A118.csv',header=TRUE,row.names=1)
genes <- read.table('A118_genes_sel200.txt',sep="\t",header=T)
#genes$Chromosome <- gsub('.[0-9]$', '', genes$Chromosome)

#samples <- c( "8577-AB +", "8577-AB -", "8580-AB +", "8580-AB -", "8581-AB +", "8581-AB -", "8582-AB +", "8582-AB -", "A118 +", "A118 -")
#samples <- sub('ipd_','',colnames(x)[grep('ipd_',colnames(x))])
samples <- colnames(x)[grep('ipd_',colnames(x))]
motifs <- unique(x$motif)

h <- rep(seq(0,0.95,by = 1.0/(perpage * 2)))
s <- rep(c(1.0,0.8),perpage * 2)
v <- rep(c(1.0,0.8),perpage * 2)
colors <- hsv(h,s,v)

l <- list(samples,c("ipd_A118","ipd_NV2_AB","ipd_NV349_AB","ipd_NV351_AB","ipd_NV353_AB","ipd_NV355_AB"),c("ipd_NV333_AB","ipd_ATCC17978"))

pdf('motifs_circos.pdf',paper="a4r",width=0,height=0)
op <- par(mar=c(1,0,1,0))

for (sl in l)
{
  # build a table of motifs counts per sample and max IPD value
  hm <- NULL
  detected <- NULL
  names <- NULL
  for (sample in sl)
  {
    xm <- NULL
    xd <- NULL
    xp <- NULL
    for(ms in motifs)
    {
      idx <- which(x$motif == ms)
      if (length(names) < length(motifs))
      {
	names <- c(names, paste(ms,x$type[idx[1]],x$modpos[idx[1]]))
      }
      idx.p <- which(x[idx,sample] != 'NA')
      mm <- 1.0
      if (length(idx.p) > 0) {
	mm <- max(x[idx.p,sample])
      }
      mt <- length(idx)
      mp <- length(idx.p)
      xm <- rbind(xm,mm)
      xd <- rbind(xd,paste(mp,'/',mt))
      xp <- rbind(xp,mp/mt)
    }
    hm <- cbind(hm,xp)
    colnames(hm)[ncol(hm)] <- sub('ipd_','',sample)
    detected <- cbind(detected,xd)
    colnames(detected)[ncol(detected)] <- sub('ipd_','',sample)
  }

  rownames(hm) <- names
  hm <- as.matrix(hm)
  rownames(detected) <- names
  detected <- as.matrix(detected)

  heatmap.2(hm,trace='n',main='Fraction motifs modified on A118',mar=c(10,15),adjCol=c(NA,0.5),sepwidth=c(0.001,0.001),sepcolor='grey20',colsep=1:ncol(hm),rowsep=1:nrow(hm),cellnote=detected,notecex=mynotecex,notecol='black',cexCol=1)

  heatmap.2(hm,trace='n',main='Fraction motifs modified on A118',mar=c(10,15),adjCol=c(NA,0.5),sepwidth=c(0.001,0.001),sepcolor='grey20',colsep=1:ncol(hm),rowsep=1:nrow(hm),cexCol=1)
}

for(motif in motifs)
{
  iAp <- which(x$motif==motif & x$strand=='+')
  iAm <- which(x$motif==motif & x$strand=='-')
  t.max <- max(x[c(iAp,iAm),samples],na.rm=T)
  t.min <- 1.0
  t.mean <- mean(as.matrix(x[c(iAp,iAm),samples]),na.rm=T)

  for (page in seq(1,length(samples),by=perpage))
  {
    cursamples <- samples[min(page,length(samples)):min(page + perpage - 1,length(samples))]
    circos.par(cell.padding = c(0.02, 0), gap.degree = c(3))
    circos.initialize(xlim = a.all)
    i <- 1
    names <- NULL
    for (sample in cursamples)
    {
      circos.track(ylim = c(t.min,t.max), track.height = 1.0/nbdiv)
      circos.trackPoints(x$chr[iAp], x=x$pos[iAp], y=x[iAp,sample],pch=1,cex=mycex,col=colors[i])
      circos.trackPoints(x$chr[iAm], x=x$pos[iAm], y=x[iAm,sample],pch=4,cex=mycex,col=colors[i+1])
      names <- c(names,paste(sub('ipd_','',sample),'+'),paste(sub('ipd_','',sample),'+'))
      i <- i + 2
    }

    legend(-1,.9,adj = c(0,0.5),xjust=0.5,yjust=0.5,cex=mycexlegend,
	   title=paste(motif,x$type[iAp[1]],x$modpos[iAp[1]], 'on A118'),
	   fill=colors[seq(1,length(names))],
	   legend=names)
    circos.yaxis(side='left',track.index=1,sector.index=rownames(a.all)[1],labels.cex=axisFontSize)

    t.gmin <- t.max + t.max / 30
    t.gmax <- t.max + t.max / 30 + t.max / 12
    glabelpos <- t.gmax + t.max / 5
    for(i in 1:nrow(a.all))
    {
      n <- rownames(a.all)[i]
      circos.axis(sector.index = n, track.index=1, major.at = seq(0,a.all[i,2],by=100000), labels.cex=axisFontSize, labels.facing='clockwise', direction='outside')
      circos.text(a.all[i,2]/2, t.max, n, track.index = 1, sector.index = n, facing = "inside", cex = axisFontSize + 0.25,adj = c(labelXpos, -labelYpos))
      for(tr in 1:length(cursamples))
      {
	circos.lines(c(1,a.all[i,2]), c(t.mean,t.mean), track.index = tr, sector.index = n, col='grey')
      }
      gi <- which(genes$Chromosome == n & genes$Strand == "+")
      if (length(gi) > 0) {
	for(tr in 1:length(cursamples))
	{
	  circos.rect(genes$Gene_Start[gi]-pad,t.gmin,genes$Gene_End[gi]+pad,t.gmax, track.index = tr, sector.index = n, col='pink', border='pink')
	}
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
		    genes$Gene_ID,
		    track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
      }
      gi <- which(genes$Chromosome == n & genes$Strand == "-")
      if (length(gi) > 0) {
	for(tr in 1:length(cursamples))
	{
	  circos.rect(genes$Gene_Start[gi]-pad,t.gmin,genes$Gene_End[gi]+pad,t.gmax, track.index = tr, sector.index = n, col='cyan', border='cyan')
	}
	circos.text(genes$Gene_Start[gi],rep(glabelpos,length(gi)),
		    genes$Gene_ID,
		    track.index=1,sector.index=n,facing='clockwise',niceFacing=TRUE,cex=axisFontSize-0.04,adj=c(0.2,0.6),col='grey30')
      }
    }

    circos.clear()
  }
}

dev.off()

#----------------------------------------------------------------------------------
