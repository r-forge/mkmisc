\name{plot.manhattan.LD}
\alias{plot.manhattan.LD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot LD pattern and MAF in a Manhattan plot.
}
\description{
Function for plotting local LD pattern (relative to a pre-selected marker) on a Manhattan plot resulting from genome-wide association study (GWAS). Each marker on the plot is colored according to its LD with the reference marker. Color codes are not continuous, but LD is discretized into LD intervals. This is useful for examining signals found in a GWAS study. In addition, minor allele frequency will be plotted in the lower panel. 
}
\usage{
plot.manhattan.LD(data, gwas.result, chr, region, index.snp, p.value = 0.05, bonferroni = T, mafThreshold = 0.05)
}
\arguments{
  \item{data}{
A \code{gwaa.data} class object as used by \code{\link[GenABEL]{gwaa.data-class}}
}
  \item{gwas.result}{
A \code{scan.gwaa} class object as used by \code{\link[GenABEL]{scan.gwaa-class}}
}
  \item{chr}{
Chromosome nubmer (name) to be displayed.
}
  \item{region}{
A vector of two coordinates to display.
}
  \item{index.snp}{
Index of the reference SNP
}
  \item{p.value}{
p-value threshold to visualize using red dashed line
}
  \item{bonferroni}{
logical indicating whether Bonferroni-correction shall be used for the displayed p-value threshold
}
  \item{mafThreshold}{
threshold value for minor allele frequency (MAF) - displayed as a red line on the lower panel.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
