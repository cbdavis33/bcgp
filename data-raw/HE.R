rm(list = ls())
cat("\014")

load("~/Documents/gitCollabos/hbcgpPaper/prediction/RData/HE.RData")

HE <- list(x = xTrain, y = as.vector(yTrain),
           xTest = xPred, yTest = as.vector(yPred))

usethis::use_data(HE, overwrite = TRUE)
