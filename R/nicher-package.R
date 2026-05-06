## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom stats AIC BIC logLik nobs
#' @useDynLib nicher, .registration = TRUE
## usethis namespace: end
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(".data")
}
