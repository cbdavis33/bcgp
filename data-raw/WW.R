rm(list = ls())
cat("\014")

load("~/Documents/gitCollabos/hbcgpPaper/prediction/RData/WW.RData")

WW <- list(x = xTrain, y = as.vector(yTrain),
           xTest = xPred, yTest = as.vector(yPred))

usethis::use_data(WW, overwrite = TRUE)
