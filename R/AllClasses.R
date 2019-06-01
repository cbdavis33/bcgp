## Contains all classes for the bcgp package
setClass(Class = "bcgpmodel",
         slots = c(data = "list",       # raw and scaled, then x and y
                   priors = "list",
                   inits = "list",
                   stationary = "logical",
                   composite = "logical",
                   algorithm = "character",
                   scaled = "character"))

setClass(Class = "bcgpfit",
         slots = c(model_pars = "character",
                   par_dims = "list",
                   sim,
                   sampler_args = "list",
                   date = "character",
                   .MISC = "environment"),
         contains = "bcgpmodel")

setClass(Class = "bcgpfitpred",
         slots = c(preds = "list"),
         contains = "bcgpfit")
