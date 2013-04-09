#' The \code{fstats.result} class
#'
#' This class contains results of running the \code{\link[cgmisc]{compute.Fstats}} function.
#'
#' This line and the next ones go into the details.
#' This line thus appears in the details as well.
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{pops}:}{A \code{list} of per-population statistics}
#'    \item{\code{glob}:}{A \code{data.frame} of global statistics}
#'  }
#'
#' @name fstats.result
#' @rdname fstats.result
#' @aliases fstats.result-class 
#' @exportClass fstats.result
#' @author Marcin Kierczak <\email{Marcin.Kierczak@@slu.se}>
setClass("fstats.result", representation(pops = "list", glob = "data.frame"))
