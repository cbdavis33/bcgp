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
#'
#'   # n_kept <- x@sim$n_save - x@sim$warmup2
#'   # cat("Inference for Stan model: ", x@model_name, '.\n', sep = '')
#'   # cat(x@sim$chains, " chains, each with iter=", x@sim$iter,
#'   #     "; warmup=", x@sim$warmup, "; thin=", x@sim$thin, "; \n",
#'   #     "post-warmup draws per chain=", n_kept[1], ", ",
#'   #     "total post-warmup draws=", sum(n_kept), ".\n\n", sep = '')
#'   #
#'   # if (!is.null(x@stan_args[[1]]$method) &&
#'   #     x@stan_args[[1]]$method == "variational") {
#'   #   print(round(s$summary, digits_summary), ...)
#'   #   cat("\nApproximate samples were drawn using VB(", x@stan_args[[1]]$algorithm, ") at ", x@date,
#'   #       ".\n", sep = '')
#'   #   message("We recommend genuine 'sampling' from the posterior distribution for final inferences!")
#'   #   return(invisible(NULL))
#'   # }
#'   #
#'   # # round n_eff to integers
#'   # s$summary[, 'n_eff'] <- round(s$summary[, 'n_eff'], 0)
#'   #
#'   # print(round(s$summary, digits_summary), ...)
#'   #
#'   # sampler <- attr(x@sim$samples[[1]], "args")$sampler_t
#'   #
#'   # cat("\nSamples were drawn using ", sampler, " at ", x@date, ".\n",
#'   #     "For each parameter, n_eff is a crude measure of effective sample size,\n",
#'   #     "and Rhat is the potential scale reduction factor on split chains (at \n",
#'   #     "convergence, Rhat=1).\n", sep = '')
#'   # return(invisible(NULL))
#' }
