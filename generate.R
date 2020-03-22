require(reshape)
require(deSolve)

require(ggplot2)


# Create a multivariate function with specified properties:
#   tmax: maximum time
#   multiplicities: number of correlations each parameter has
#   degrees: polynomial degree of each parameter
#   dimension: the dimension of the output
#   returns a multivariate function of the vector of parameters and time
makeGenerator <- function(multiplicities, degrees, dimension) {

    single <- function(degree) {
      x0 <- runif(1)
      z0 <- runif(1)
      function(x) {
          if (x < x0)
              0
          else
              z0 * (x - x0)^degree
      }
    }
  
    locations <- lapply(multiplicities, function(m) sample(1:dimension, m))
    functions <- lapply(degrees, single)
    
    start <- runif(dimension, -0.25, 0.75)
    coefs <- matrix(runif(dimension^2, -0.25, 0.75), dimension, dimension)    
    shift <- matrix(runif(dimension^2, -0.25, 0.75), dimension, dimension)
    
    function(x, ts) {
        z <- rep(0, dimension)
        for (i in 1:length(locations))
            for (j in locations[[i]])
                z[j] <- z[j] + functions[[i]](x[i])
        ode(start, ts, function(t, y, params) {list((coefs %*% y) * z * (1 - ((shift %*% y) * z)))})
    }
    
}


### Example

# Create an example function with three input dimensions and three output dimensions.
# The first two input parameters affect two outputs, and the third affects all outputs.
# The first parameter involves a discontinuity in the output, the second parameter involves a discontinuity in the derivative,
# and the third parameter involves a discontinuity in the second derivative.

f <- makeGenerator(c(2, 2, 3), c(0, 1, 2), 3)

# Evaluate at some times.

ts <- seq(0, 10, 0.5)

# Run a single simulation.

runCase <- function(case) {
  x <- runif(3)
  y <- f(x, ts)
  data.frame(case=factor(case), x1=x[1], x2=x[2], x3=y[3], t=ts, y1=y[, 2], y2=y[, 3], y3=y[, 4])
}

# Run multiple cases.

ys <- NULL
for (i in 1:10)
  ys <- rbind(ys, runCase(i))

# Make some plots.

ggplot(ys, aes(x=t, y=y1, color=case)) + geom_line()

ggplot(ys, aes(x=t, y=y2, color=case)) + geom_line()

ggplot(ys, aes(x=t, y=y3, color=case)) + geom_line()
