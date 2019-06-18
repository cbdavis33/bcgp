setMethod("initialize", "bcgpsims",
          function(.Object, ...){
              data <- simulate_from_model()
              toReturn@data <- with(data@data, list(x = x, y = y))
              toReturn@pred <- with(data@data, list(x = xPred, y = yPred))
              parameters <- data@parameters
              stationary <- data@stationary
              composite <- data@composite
              seed <- data@seed
              return(toReturn)
            })

# setMethod("plot", signature = "bcgpsims",
#           function(object){
#
#           })
