rm(list = ls())
cat("\014")

load("~/Documents/gitCollabos/hbcgpPaper/prediction/RData/HE.RData")

HE <- list(x = xTrain, y = as.vector(yTrain),
           xPred = xPred, yPred = as.vector(yPred))

usethis::use_data(HE, overwrite = TRUE)
