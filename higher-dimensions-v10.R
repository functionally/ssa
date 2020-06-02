
require(data.table)
require(deSolve)
require(magrittr)
require(randtoolbox)

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
  
    locations <- lapply(multiplicities, function(m) sample(1:dimension, m, replace=TRUE))
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


RNGkind("Mersenne-Twister", "Inversion", "Rejection")

set.seed(46)

n <- 10000

ts <- seq(0, 3, 0.2)
ts

for (m in c(3, 20))
    for (k in c(3, 20))
        if (m != 3 || k != 3) {
            f.multiplicities <- sample(1:5, m, replace=TRUE)
            f.degrees <- sample(0:3, m, replace=TRUE)
            f <- makeGenerator(f.multiplicities, f.degrees, k)
            xs <- sobol(n, m)
            colnames(xs) <- paste("x", 1:m, sep="")
            rownames(xs) <- 1:n
            ys <- NULL
            for (i in 1:n) {
                x <- xs[i, ]
                if (m == 3)
                    y <- f(x, 10 * ts)
                else
                    y <- f(x, ts)
                ys <- rbind(ys, data.table(case=i, y))
            }
            xs <- data.table(case = 1:n, xs)
            colnames(ys) <- c("case", "time", paste("y", 1:k, sep=""))
            write.table(xs, file=paste("xs-", m, "x", k, "-v10.csv", sep=""), sep=",", quote=FALSE, row.names=FALSE)
            write.table(ys, file=paste("ys-", m, "x", k, "-v10.csv", sep=""), sep=",", quote=FALSE, row.names=FALSE)
        }

fread("ys-3x20-v10.csv") %>% summary

fread("ys-20x3-v10.csv") %>% summary

fread("ys-20x20-v10.csv") %>% summary
