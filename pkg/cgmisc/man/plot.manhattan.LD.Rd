\name{plot.manhattan.LD}
\alias{plot.manhattan.LD}
\title{Plot LD pattern and MAF in a Manhattan plot}
\usage{
  plot.manhattan.LD(data, gwas.result, chr, region,
    index.snp, p.value = 0.05, bonferroni = T,
    mafThreshold = 0.05)
}
\arguments{
  \item{data}{a gwaa.data class object as used by
  \code{\link[GenABEL]{gwaa.data-class}}}

  \item{gwas.result}{a scan.gwaa class object as used by
  \code{\link[GenABEL]{scan.gwaa-class}}}

  \item{chr}{chromosome nubmer (name) to be displayed}

  \item{region}{a vector of two coordinates to display}

  \item{index.snp}{index of the reference SNP}

  \item{p.value}{p-value threshold to visualize using red
  dashed line}

  \item{bonferroni}{logical indicating whether
  Bonferroni-correction shall be used for the displayed
  p-value threshold}

  \item{mafThreshold}{threshold value for minor allele
  frequency (MAF) -- displayed as a red line on the lower
  panel}
}
\value{
  NULL
}
\description{
  Function for plotting local LD pattern (relative to a
  pre-selected marker) on a Manhattan plot resulting from
  genome-wide association study (GWAS). Each marker on the
  plot is colored according to its LD with the reference
  marker. Color codes are not continuous, but LD is
  discretized into LD intervals. This is useful for
  examining signals found in a GWAS study. In addition,
  minor allele frequency will be plotted in the lower
  panel.
}
\examples{
\dontrun{
plot.manhattan.LD()
}
}
\author{
  Marcin Kierczak <\email{Marcin.Kierczak@slu.se}>
}
\seealso{
  \code{\link[GenABEL]{gwaa.data-class}},
  \code{\link[GenABEL]{scan.gwaa-class}}
}
\keyword{LD}
\keyword{Manhattan}
\keyword{disequilibrium}
\keyword{genetics}
\keyword{linkage}

