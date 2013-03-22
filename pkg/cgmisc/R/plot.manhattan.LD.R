plot.manhattan.LD <-
function(data, gwas.result, chr, region, index.snp, p.value=0.05, bonferroni=T, mafThreshold=.05) {
  shift <- 2.3
  topMargin <- 0
  # Helper function for retrieving MAF in a region (chromosome)
  getMAF <- function(data, region) {
    summ <- summary(gtdata(data)[,region])
    tmp_maf <- (summ$P.22 + 0.5 * summ$P.12)/summ$NoMeasured
    tmp_ma_count <- summ$P.22 + 0.5 * summ$P.12
    tmp_het <- summ$P.12 / summ$NoMeasured
    maf <- data.frame(snp=as.character(row.names(summ)), 
                      chr=summ$Chromosome, pos=map(gtdata(data)[,region]), 
                      allele=summ$A2, maf=tmp_maf, maCnt=tmp_ma_count, heterozygosity=tmp_het)
  }
  startCoord <- region[1]
  stopCoord  <- region[2]
  myChromosome <- data@gtdata[ ,which(data@gtdata@chromosome == chr)]
  region <- which(myChromosome@map >= startCoord & myChromosome@map <= stopCoord)
  r2matrix <- r2fast(myChromosome, snpsubset = region)
  r2matrix[lower.tri(r2matrix)] <- t(r2matrix)[lower.tri(r2matrix)]
  markers <- which(data@gtdata@snpnames %in% names(region)) 
  markers.coords <- data@gtdata@map[markers]
  idx.marker <- which(data@gtdata@snpnames == index.snp)
  idx.marker.coords <- data@gtdata@map[idx.marker]
  pvals <- -log10(gwas.result@results$P1df[markers])
  plot(markers.coords, pvals, type='n', xlab="Position (Mb)", ylab=expression(-log[10](p-value)), ylim=c(-shift, max(pvals)+3), axes=F)
  # Plot grid
  abline(h=seq(.5, max(pvals) + topMargin, .5), col="grey", lty=3)
  #points(idx.marker.coords, -log10(gwas.mm@results$P1df[idx.marker]), col="red", pch=19, cex=1)
  r2vec <- r2matrix[index.snp,]
  r2vec[is.na(r2vec)] <- -1
  r2col <- cut(r2vec, breaks=c(1.0,0.8,0.6,0.4,0.2,0.0,-1), 
               labels=rev(c("#9E0508","tomato","chartreuse3","cyan3","navy", "black")), include.lowest=T)
  r2pch <- rep(19, length(r2col))
  r2pch[which(r2col == "black")] <- 1
  points(markers.coords, -log10(gwas.result@results$P1df[markers]), col=as.character(r2col), pch=r2pch, cex=.8)
  if (bonferroni) {
    p.value <- -log10(p.value/nids(data))
    abline(h=-log10(p.value), col="red", lty=2)
  }
  legend(startCoord + 10, max(pvals) + 1, legend=c("0.8 - 1.0","0.6 - ", "0.4 -", "0.2 -", "0.0 -"), pch=19, bty='n', 
         col=c("#9E0508","tomato","chartreuse3","cyan3","navy"), cex=.7, title=expression(r^2))
  # Plot MAF below Manhattan (shift units lower)
  maf <- getMAF(data, region)
  lines(markers.coords, 4 * maf$maf - shift, col="#2C7FB8")
  # Plot horizontal line at MAF = 0.05
  abline(h=mafThreshold * 4 - shift, col=rgb(1,0,0,1), lty=2)
  # Plot axes
  step <- (stopCoord - startCoord) / 5
  axis(1, at = seq(startCoord, stopCoord, by=step), labels=format(seq(startCoord, stopCoord, by=step)/1e6, scientific=F, digits=3))
  axis(2, at = 0:(max(pvals) + topMargin + 1))
  axis(4, at = c(2 - shift, 1 - shift, 0.2 - shift, 0 - shift), labels = c(.5, .25, .05, 0))
  mtext("MAF", side=2, at=0.2-shift, outer=F)
  
}
