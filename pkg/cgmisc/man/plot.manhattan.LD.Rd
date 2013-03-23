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
}
\references{
}

\author{
Marcin Kierczak
}

\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mnhattan plot }
\keyword{ GWAS }
\keyword{ LD }
\keyword{ linkage disequilibrium}