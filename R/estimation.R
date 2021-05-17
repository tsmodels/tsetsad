#' Estimates an ETS model given a specification object using maximum likelihood and autodiff
#'
#' @param object An object of class tsets.spec.
#' @param solver One of either \dQuote{nloptr} or \dQuote{nlminb}. Using nlminb which takes
#' the hessian appears to provide better results.
#' @param control Solver control parameters.
#' @param ... additional parameters passed to the estimation function
#' @return An list of coefficients and other information.
#' @details This function is not expected to be used by itself but rather as a plugin
#' to be called from the estimate method of the tsets package.
#' @export estimate_ad.tsets.spec
#' @aliases estimate_ad
#' @export
estimate_ad.tsets.spec <- function(object, solver = "nlminb", control = list(trace = 0), ...)
{
    fun <- switch(as.character(object$model$type),
                  "1" = estimate_aaa_ad,
                  "2" = estimate_mmm_ad,
                  "3" = estimate_mam_ad,
                  "4" = estimate_powermam_ad)
    solution <- fun(object, solver = solver, control = control, ...)
    return(solution)
}