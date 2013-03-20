\name{plot.manhattan.LD}
\alias{plot.manhattan.LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot LD and MAF on Manhattan plot
}
\description{
This is a function for ploting color-coded linkage disequilibrium on a Manhattan plot. LD is relative to a selected marker. In addition, minor allele frequency (MAF) will be plotted below Manhattan. 
}
\usage{
plot.manhattan.LD(data, gwas.result, chr, region, index.snp, p.value = 0.05, bonferroni = T, mafThreshold = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an object of the gwaa.data class (GenABEL internal class)
}
  \item{gwas.result}{
an object of the gwaa.scan class containing association testing results (GenABEL class)
}
  \item{chr}{
the chromosome to be plotted
}
  \item{region}{
a vector of genomic coordinates to be plotted (chromosome-wise coords)
}
  \item{index.snp}{
a reference marker. LD will be computed with respect to this marker.
}
  \item{p.value}{
a p-value threshold for plotting a threshold line
}
  \item{bonferroni}{
a logical indicating whether Bonferroni correction shall be used
}
  \item{mafThreshold}{
a threshold for plotting MAF threshold line
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
no return value, a plot will be displayed
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Marcin Kierczak
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, gwas.result, chr, region, index.snp, p.value = 0.05, 
    bonferroni = T, mafThreshold = 0.05) 
{
    shift <- 2.3
    topMargin <- 0
    getMAF <- function(data, region) {
        summ <- summary(gtdata(data)[, region])
        tmp_maf <- (summ$P.22 + 0.5 * summ$P.12)/summ$NoMeasured
        tmp_ma_count <- summ$P.22 + 0.5 * summ$P.12
        tmp_het <- summ$P.12/summ$NoMeasured
        maf <- data.frame(snp = as.character(row.names(summ)), 
            chr = summ$Chromosome, pos = map(gtdata(data)[, region]), 
            allele = summ$A2, maf = tmp_maf, maCnt = tmp_ma_count, 
            heterozygosity = tmp_het)
    }
    startCoord <- region[1]
    stopCoord <- region[2]
    myChromosome <- data@gtdata[, which(data@gtdata@chromosome == 
        chr)]
    region <- which(myChromosome@map >= startCoord & myChromosome@map <= 
        stopCoord)
    r2matrix <- r2fast(myChromosome, snpsubset = region)
    r2matrix[lower.tri(r2matrix)] <- t(r2matrix)[lower.tri(r2matrix)]
    markers <- which(data@gtdata@snpnames \%in\% names(region))
    markers.coords <- data@gtdata@map[markers]
    idx.marker <- which(data@gtdata@snpnames == index.snp)
    idx.marker.coords <- data@gtdata@map[idx.marker]
    pvals <- -log10(gwas.result@results$P1df[markers])
    plot(markers.coords, pvals, type = "n", xlab = "Position (Mb)", 
        ylab = expression(-log[10](p - value)), ylim = c(-shift, 
            max(pvals) + 3), axes = F)
    abline(h = seq(0.5, max(pvals) + topMargin, 0.5), col = "grey", 
        lty = 3)
    r2vec <- r2matrix[index.snp, ]
    r2vec[is.na(r2vec)] <- -1
    r2col <- cut(r2vec, breaks = c(1, 0.8, 0.6, 0.4, 0.2, 0, 
        -1), labels = rev(c("#9E0508", "tomato", "chartreuse3", 
        "cyan3", "navy", "black")), include.lowest = T)
    r2pch <- rep(19, length(r2col))
    r2pch[which(r2col == "black")] <- 1
    points(markers.coords, -log10(gwas.result@results$P1df[markers]), 
        col = as.character(r2col), pch = r2pch, cex = 0.8)
    if (bonferroni) {
        p.value <- -log10(p.value/nids(data))
        abline(h = -log10(p.value), col = "red", lty = 2)
    }
    legend(startCoord + 10, max(pvals) + 1, legend = c("0.8 - 1.0", 
        "0.6 - ", "0.4 -", "0.2 -", "0.0 -"), pch = 19, bty = "n", 
        col = c("#9E0508", "tomato", "chartreuse3", "cyan3", 
            "navy"), cex = 0.7, title = expression(r^2))
    maf <- getMAF(data, region)
    lines(markers.coords, 4 * maf$maf - shift, col = "#2C7FB8")
    abline(h = mafThreshold * 4 - shift, col = rgb(1, 0, 0, 1), 
        lty = 2)
    step <- (stopCoord - startCoord)/5
    axis(1, at = seq(startCoord, stopCoord, by = step), labels = format(seq(startCoord, 
        stopCoord, by = step)/1e+06, scientific = F, digits = 3))
    axis(2, at = 0:(max(pvals) + topMargin + 1))
    axis(4, at = c(2 - shift, 1 - shift, 0.2 - shift, 0 - shift), 
        labels = c(0.5, 0.25, 0.05, 0))
    mtext("MAF", side = 2, at = 0.2 - shift, outer = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Manhattan plot }
\keyword{ LD }
\keyword{ MAF }
\keyword{ GWAS }