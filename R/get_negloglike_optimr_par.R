#' get_negloglike_optim_par
#' Function  to optime negloglike function given parameters from presence points and a sample
#' of environmental points from a M hypothesis. Internally runs Returns a list of parameters
#'
#' @param M_pts A dataframe with a sample of environmental values inside an M hypothesis region.
#' @param env_pts A dataframe with environmental variables extracted from presence points
#' @param lower logical if is set create a low bowndarie for parameters from M points
#' @param itnmax The maximum number of iterations.
#' @param fastmode logical if is true the L-BFGS-B method is by default
#' @param return_best logical. If true return the best optimization. If false return all the possible optimizations
#' @importFrom  utils as.relistable relist
#' @importFrom optimx optimx
#' @return A list with optimized parameters for negloglike function
#' @export
#'
#' @examples
#' get_negloglike_optimr_par(head(spOccPnts, 10), samMPts, fastmode = TRUE, itnmax = 1)
get_negloglike_optimr_par <- function(env_pts, M_pts, lower = FALSE, itnmax = 100, fastmode = FALSE, return_best = TRUE){
  par <- get_ellip_par(env_pts)
  initial.param <- as.relistable(par)
  ul <- unlist(initial.param)
  guess <- relist(ul)
  like.fn <- function(param.vector)
  {
    param <- relist(param.vector, skeleton = par)
    negloglike_multivariable(param$mu, param$S, sam1 = env_pts, sam2 = M_pts)
  }

  suppressWarnings({


  if(lower == TRUE){
    cat("fastmode\n")
    cat(fastmode)
    cat("\n")

    if(fastmode == FALSE){
      methods <- c("L-BFGS-B", "bobyqa", "nlminb", "nlm", "Rcgmin", "Rvmmin")
    } else {
      methods <- c("L-BFGS-B")
    }


    min.par <- get_ellip_par(M_pts)

    cat("Optimization starting \n")
    cat("Methods used:\n")
    cat(
      paste0(methods, collapse = " \n")
    )
    cat("\n")


    find.mle <- optimx::optimx(par = unlist(par),
                                 fn = like.fn,
                                 lower = unlist(par),
                                 method = methods,
                                 itnmax = itnmax)
    #})


  } else {
    cat("fastmode\n")
    cat(fastmode)
    cat("\n")



    if(fastmode == FALSE){
      methods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")
    } else {
      methods <- c("L-BFGS-B")
    }

    cat("Optimization starting \n")
    cat("Methods used:\n")
    cat(
      paste0(methods, collapse = " \n")
    )
    cat("\n")


    find.mle <- optimx::optimx(par = unlist(par),
                                 fn = like.fn,
                                 method= methods,
                                 itnmax = itnmax)
    }
  })

  mle.par <- Map(function(x)  relist(unlist(find.mle[x, 1:6]), skeleton = par), 1:nrow(find.mle))
  names(mle.par) <- rownames(find.mle)

  if(return_best == TRUE){
    simVector <- Reduce(c, Map(function(x) isSymmetric(round(x$S, 5) ), mle.par))

    find.mle <- find.mle[simVector, ]
    find.mle <- find.mle[find.mle$value == min(find.mle$value), ]

    optim_parameters <- unlist(find.mle[1, 1:6])
    optim_parameters <- round(optim_parameters, 5)
    optim_parameters <- relist(optim_parameters, skeleton = par)

    cat("Optimization finished. \n")
    cat(rownames(find.mle))
    cat(" is the returned method")
    cat("\n")
    return(optim_parameters)

  } else {
    cat("Optimization finished. \n")
    cat("All methods are returned")
    return(mle.par)
  }
}
