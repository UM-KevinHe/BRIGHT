#' GWAS to correlation
#'
#' Recover marginal correlation information from GWAS summary statistics.
#'
#' @param p a numeric vector of the p values from marginal SLR.
#' @param n a numeric vector of the sample sizes.
#' @param sign a numeric vector of the signs of coefficient estimates.
#' @return Recovered XtY/n from the GWAS summary.
#' @examples 
#' p = c(0.01,0.015)
#' n = c(100,1000)
#' sign = c(-1,1)
#' p2cor(p, n, sign)
p2cor <- function(p, n, sign = rep(1, length(p))) {
  
  stopifnot(length(n)==1 || length(n) == length(p))
  stopifnot(length(p) == length(sign))
  
  t <- sign(sign) * qt(p/2, df = n-2, lower.tail = F)
  
  return(t / sqrt(n - 2 + t^2))
  
}