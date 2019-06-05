## Contains all classes for the bcgp package
setClass(Class = "bcgppriors",
         slots = c(priors = "list",
                   stationary = "logical",
                   composite = "logical",
                   noise = "logical"))

setClass(Class = "bcgpinits",
         slots = c(inits = "list"),
         contains = "bcgppriors")

setClass(Class = "bcgpmodel",
         slots = c(data = "list",       # raw and scaled, then x and y
                   priors = "list",
                   inits = "list",
                   stationary = "logical",
                   composite = "logical",
                   algorithm = "character",
                   scaled = "logical"))

setClass(Class = "bcgpfit",
         slots = c(model_pars = "character",
                   par_dims = "list",
                   sim = "list",
                   sampler_args = "list",
                   date = "character",
                   .MISC = "environment"),
         contains = "bcgpmodel")

setClass(Class = "bcgpfitpred",
         slots = c(preds = "list"),
         contains = "bcgpfit")
