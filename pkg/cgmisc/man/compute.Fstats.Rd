\name{compute.Fstats}
\alias{compute.Fstats}
\title{
Function for computing \eqn{F_{st}}{F[st]} fixation index for two populations.
}
\description{
Fixation index (\eqn{F_{st}}{F[st]}) is a measure of population differentiation due to genetic structure.
Given a set of genotypes in two populations, the function computes fixation index Fst and corresponding indices: \eqn{F_{it}}{F[it]} and \eqn{F_{is}}{F[is]}. 
}
\usage{
compute.Fstats(data, pops)
}

\arguments{
	\itemize{
  	\item data -- a \code{\link[GenABEL]{gwaa.data-class}} object,
  	\item pops -- a vector of two values indicating to which population an individual belongs.	
	}
}

\details{
	\itemize{
  	\item data -- a standard \code{\link[GenABEL]{gwaa.data-class}} object.
  	\item pops -- a vector of two values indicating to which population an individual belongs. Typically, one uses a vector of zeroes and ones where 0 marks an individual belonging to population 1 and 1 marks an individual belonging to population 2. Often, the vector is a result of clustering in MDS-scaled genomic kinship space. 	
	}
}

\value{
	An S4 class object with two slots:
	\itemize{
		\item \code{@pops} -- a list of two data frames containing per population allele frequency statistics,
		\item \code{@glob} -- a data frame containing global statistics.
	}
	Each \code{@pops} element contains the following columns per every marker (row):
	\itemize{
		\item N -- number of non-missing genotypes,
		\item p -- major allele frequency,
		\item q -- minor allele frequency,
		\item exp.AA -- expected count of major allele homozygotes,
		\item exp.Aa -- expected count of heterozygotes,
		\item exp.aa -- expected count of minor allele homozygotes,
		\item obs.het -- observed heterozygosity,
		\item exp.het -- expected heterozygosity,
		\item dev.het -- deviation from expected heterozygosity,
		\item Fs - inbreeding coefficient.
	}

	The \code{@glob} data frame contains global (joint) per marker information in columns:
	\itemize{
		\item p -- global major allele frequency,
		\item q -- global minor allele frequency,
		\item HI -- global heterozygosity index HI,
		\item HS -- global heterozygosity index HS,
		\item HT -- global heterozygosity index HT,
		\item FIS -- \eqn{F_{is}}{F[is]} statistics,
		\item FST -- \eqn{F_{st}}{F[st]} statistics,
		\item FIT -- \eqn{F_{fit}}{F[fit]} statistics.
	}
}
\references{
Holsinger, Kent E.; Bruce S. Weir (2009). "Genetics in geographically structured populations: defining, estimating and interpreting FST". Nat Rev Genet 10 (9): 639–650. doi:10.1038/nrg2611. ISSN 1471-0056. PMID 19687804.
}
\author{
Marcin Kierczak
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[GenABEL]{gwaa.data-class}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FST }
\keyword{ fixation index }
\keyword{ population structure }
\keyword{ heterozygosity }
