\name{clump.markers}
\alias{clump.markers}
\title{Clump markers according to their LD}
\usage{
  clump.markers(data, an, chr = 1, bp.dist = 250000,
    p1 = 1e-04, p2 = 0.01, r2 = 0.5, image = F)
}
\arguments{
  \item{p1}{threshold for index markers,}

  \item{p2}{threshold for clumping,}

  \item{r2}{threshold for LD,}

  \item{bp.dist}{threshold for inter-marker distance,}

  \item{chr}{chromosome to be clumped,}

  \item{an}{\code{\link[GenABEL]{gwaa.scan-class}}
  objectwith association test results,}

  \item{data}{data object in
  \code{\link[GenABEL]{gwaa.data-class}},}

  \item{image}{a logical indicating whether to plot
  clumping results or not.}
}
\value{
  a list of clumps
}
\description{
  \code{clumpMarkers} implements clumping procedure
  described in PLINK documentation on a
  \code{\link[GenABEL]{gwaa.data-class}} object.
}
\examples{
\dontrun{
   clumps <- clumpMarkers(data.qc0, image=T, an <- an0, chr <- 6, bp.dist <- 250e3, p1 <- 0.0001, p2 <- 0.01, r2 <- 0.5)
   }
}
\author{
  Marcin Kierczak <\email{Marcin.Kierczak@slu.se}>
}
\references{
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/clump.shtml}
}

