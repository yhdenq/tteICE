#' @title
#' Plot method for 'tteICE' objects
#'
#' @description
#' This function plots the estimated potential cumulative incidence functions or treatment effect curve
#' with pointwise confidence intervals.
#'
#' @param x
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @param type
#' Which plot to create: \code{type="ate"} indicates to plot the estimated treatment effects;
#' \code{type="inc"} indicates to plot the estimated cumulative incidence functions (CIFs).
#'
#' @param decrease
#' Corresponds to the argument in \code{\link[tteICE]{plot_ate}} and \code{\link[tteICE]{plot_inc}}.
#'
#' @param conf.int
#' #' Confidence level for the pointwise confidence intervals
#' If \code{conf.int = NULL}, no confidence intervals are provided.
#'
#' @param xlab
#' Label for the x-axis.
#'
#' @param xlim
#' A numeric vector of length 2 specifying the limits of the x-axis.
#' If \code{xlim=NULL} (default), the range is determined automatically from the data.
#'
#' @param ylim
#' A numeric vector of length 2 giving the limits of the y-axis.
#' If \code{ylim=NULL} (default), the range is determined automatically by the type of plot,
#' corresponding to the argument in \code{\link[tteICE]{plot_ate}} and \code{\link[tteICE]{plot_inc}}.
#'
#' @param plot.configs
#' A named \code{list} of additional plot configurations. See details in \code{\link{plot_ate}} and \code{\link{plot_inc}}
#'
#' @param ... Other arguments in function \code{\link{plot.default}} or function \code{\link{curve}}
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' bmt$A = A
#'
#' ## simple model fitting and plotting
#' fit1 = tteICE(Surv(t2, factor(d4))~A, data=bmt)
#' plot(fit1, type="ate")
#' plot(fit1, type="inc")
#'
#'
#' ## plot cumulative incidence functions with p-values
#' fit2 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
#' plot(fit2, type="inc", decrease=TRUE, ylim=c(0,1),
#'      plot.configs=list(show.p.value=TRUE))
#'
#' ## plot treatment effects for semicompeting risk data
#' fit3 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#' plot(fit3, type="ate", ylim=c(-1,1), xlab="time",
#'      plot.configs=list(col="red"))
#'
#' @seealso
#' \code{\link[tteICE]{plot_ate}}, \code{\link[tteICE]{plot_inc}},
#' \code{\link[tteICE]{surv.tteICE}}, \code{\link[tteICE]{scr.tteICE}},
#' \code{\link[tteICE]{tteICE}}
#'
#'
#' @method plot tteICE
#' @return Plot the results from a tteICE object
#' @export

plot.tteICE <- function(x, type=c("ate","inc")[1],
  decrease=FALSE,conf.int=.95,xlab='Time',xlim=NULL, ylim=NULL,
  plot.configs=list(),...){

  # if (!inherits(x, "tteICE"))  stop("`fit` must be an object returned by `surv.tteICE` or `scr.tteICE`.", call. = FALSE)
  type <- match.arg(type, c("ate", "inc"))

  if(type=="ate") {
    if(is.null(ylim)) ylim=c(-1,1) else ylim=ylim
    plot_ate(fit=x,decrease=decrease,conf.int=conf.int,xlab=xlab,xlim=xlim,ylim=ylim,
      plot.configs=plot.configs,...)
    }
  else {
    if(is.null(ylim)) ylim=c(0,1) else ylim=ylim
    plot_inc(fit=x,decrease=decrease,conf.int=conf.int,xlab=xlab, xlim=xlim,ylim=ylim,
      plot.configs=plot.configs,...)
    }
}
