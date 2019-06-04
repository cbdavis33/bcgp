rm(list = ls())
cat("\014")

load("~/Documents/gitCollabos/hbcgpPaper/prediction/RData/WW.RData")

WW <- list(x = xTrain, y = as.vector(yTrain),
           xPred = xPred, yPred = as.vector(yPred))

usethis::use_data(WW, overwrite = TRUE)
