#' #' @export
#' setMethod("show", signature(x = "bcgpfit"),
#'           function(x){
#'             print.bcgpfit(x, pars = x@model_pars)
#'
#'           })
#'
#' print.bcgpfit <- function(x, pars = x@model_pars,
#'               probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
#'               digits_summary = 2){
#'
#'   s <- summary(x, pars, probs)
#'   if (is.null(s)) return(invisible(NULL))
#' }
