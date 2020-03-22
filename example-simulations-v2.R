
require(data.table)
require(deSolve)
require(magrittr)

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


runCases <- function(f, xs, ts) {
    result <- NULL
    for (i in 1:nrow(xs))
        result <- rbind(
            result,
            data.table(
                case=xs[i, case],
                f(as.numeric(xs[i, .(x1, x2, x3)]), ts)
            )
        )
    result[, .(
        case,
        t=time,
        y1=`1`,
        y2=`2`,
        y3=`3`
    )]
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

set.seed(46)

f <- makeGenerator(c(2, 2, 3), c(0, 1, 2), 3)

ts <- seq(0, 10, 0.5)

xs <- data.table(
        cbind(
        case=1:11^3,
        expand.grid(
            x1=seq(0, 1, 0.1),
            x2=seq(0, 1, 0.1),
            x3=seq(0, 1, 0.1)
        )
    )
)
xs %>% dim

xs %>% summary

ys <- runCases(f, xs, ts)
ys %>% dim

ys %>% summary

ggplot(ys, aes(x=t, y=y1, group=case)) + geom_line()

ggplot(ys, aes(x=t, y=y2, group=case)) + geom_line()

ggplot(ys, aes(x=t, y=y3, group=case)) + geom_line()

write.csv(xs, file="xs-3d-20200322a.csv", row.names=FALSE)
write.csv(ys, file="ys-3d-20200322a.csv", row.names=FALSE)

xs2d <- xs[x3 == 0.8, .(case, x1, x2)]
xs2d %>% dim

ys2d <- merge(xs2d, ys, on=case)[, .(case, t, y1, y2)]

ys2d %>% dim

ggplot(ys2d, aes(x=t, y=y1, group=case)) + geom_line()

ggplot(ys2d, aes(x=t, y=y2, group=case)) + geom_line()

write.csv(xs2d, file="xs-2d-20200322a.csv", row.names=FALSE)
write.csv(ys2d, file="ys-2d-20200322a.csv", row.names=FALSE)
