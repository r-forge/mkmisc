\name{compute.Fstats}
\alias{compute.Fstats}
\title{
Function for computing Fst fixation index for two populations.
}
\description{
Fixation index (Fst) is a measure of population differentiation due to genetic structure.
Given a set of genotypes in two populations, the function computes fixation index Fst and corresponding indices: Fit and Fis. 
}
\usage{
compute.Fstats(data, pops)
}

\arguments{
  \item{data}{
	a \code{\link[GenABEL]{gwaa.data-class} object
}
  \item{pops}{
	vector of two values indicating to which population an individual belongs	
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An S4 class object with two slots:
\item{@pops }{a list of two data frames containing per population allele frequency statistics.}
\item{@glob }{a data frame containing global statistics. }

Each \code{@pops} element contains the following columns per every marker (row):
\item{N}{number of non-missing genotypes}
\item{p}{major allele frequency}
\item{q}{minor allele frequency}
\item{exp.AA}{expected count of major allele homozygotes}
\item{exp.Aa}{expected count of heterozygotes}
\item{exp.aa}{expected count of minor allele homozygotes}
\item{obs.het}{observed heterozygosity}
\item{exp.het}{expected heterozygosity}
\item{dev.het}{deviation from expected heterozygosity}
\item{Fs}{inbreeding coefficient}

The \code{@glob} data frame contains global (joint) per marker information in columns:
\item{p}{global major allele frequency}
\item{q}{global minor allele frequency}
\item{HI}{global heterozygosity index HI}
\item{HS}{global heterozygosity index HS}
\item{HT}{global heterozygosity index HT}
\item{FIS}{F_is statistics}
\item{FST}{F_st statistics}
\item{FIT}{F_it statistics}
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FST }
\keyword{ fixation index }
\keyword{ population structure }
\keyword{ heterozygosity }
