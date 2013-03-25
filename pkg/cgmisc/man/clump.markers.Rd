\name{clump.markers}
\alias{clump.markers}
\title{Clump markers according to their LD
\code{clumpMarkers} implements clumping procedure described in PLINK documentation on
a \code{\link[GenABEL]{gwaa.data-class}} object.}
\usage{
  clump.markers(data, an, chr = 1, bp.dist = 250000,
    p1 = 1e-04, p2 = 0.01, r2 = 0.5, image = F)
}
\arguments{
  \item{p1}{threshold for index markers}

  \item{p2}{threshold for clumping}

  \item{r2}{threshold for LD}

  \item{bp.dist}{threshold for inter-marker distance}

  \item{chr}{chromosome to be clumped}

  \item{an}{a \code{\link[GenABEL]{gwaa.scan-class}} object
  with association test results}

  \item{data}{data object in
  \code{\link[GenABEL]{gwaa.data-class}}}

  \item{image}{a logical indicating whether to plot
  clumping results or not}
}
\value{
  a list of clumps
}
\description{
  Clump markers according to their LD \code{clumpMarkers}
  implements clumping procedure described in PLINK
  documentation on a \code{\link[GenABEL]{gwaa.data-class}}
  object.
}
\author{
  Marcin Kierczak \email{Marcin.Kierczak
}
\references{
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/clump.shtml}
}

