
require(abind)
require(data.table)
require(magrittr)
require(rTensor)

require(ggplot2)
require(GGally)

xs3d <- fread("xs-3d-20200322a.csv")[, case:=factor(case)]
xs3d %>% dim

ys3d <- fread("ys-3d-20200322a.csv")[, case:=factor(case)]
ys3d %>% dim

s3d <- merge(xs3d, ys3d)
s3d %>% dim

y <- array(
    as.matrix(s3d[order(t, x3, x2, x1)][, .(y1, y2, y3)]),
    dim = c(
        length(unique(s3d$x1)),
        length(unique(s3d$x2)),
        length(unique(s3d$x3)),
        length(unique(s3d$t )),
        3
    ),
    dimnames = list(
        x1=sort(unique(s3d$x1)),
        x2=sort(unique(s3d$x2)),
        x3=sort(unique(s3d$x3)),
        t =sort(unique(s3d$t )),
        i =1:3
    )
)
y %>% dim

dy1 <- y[2:11,     ,     , , ] - y[1:10,     ,     , , ]
dy2 <- y[    , 2:11,     , , ] - y[    , 1:10,     , , ]
dy3 <- y[    ,     , 2:11, , ] - y[    ,     , 1:10, , ]

dy <- abind(
    (dy1[    , 2:11, 2:11, , ] + dy1[    , 1:10, 1:10, , ]) / 2,
    (dy2[2:11,     , 2:11, , ] + dy2[1:10,     , 1:10, , ]) / 2,
    (dy3[2:11, 2:11,     , , ] + dy3[1:10, 1:10,     , , ]) / 2,
    along=6
)
dy %>% dim

ddy1 <- dy[2:10,     ,     , , , ] - dy[1:9,    ,    , , , ]
ddy2 <- dy[    , 2:10,     , , , ] - dy[   , 1:9,    , , , ]
ddy3 <- dy[    ,     , 2:10, , , ] - dy[   ,    , 1:9, , , ]

ddy <- abind(
    (ddy1[    , 2:10, 2:10, , , ] + ddy1[   , 1:9, 1:9, , , ]) / 2,
    (ddy2[2:10,     , 2:10, , , ] + ddy2[1:9,    , 1:9, , , ]) / 2,
    (ddy3[2:10, 2:10,     , , , ] + ddy3[1:9, 1:9,    , , , ]) / 2,
    along=7
)
ddy %>% dim

cy <- array(
    0,
    dim = c(
        dim(ddy)[1],
        dim(ddy)[2],
        dim(ddy)[3],
        dim(ddy)[4],
        dim(ddy)[5],
        3
    ),
    dimnames = list(
        x1=dimnames(y)$x1[2:10],
        x2=dimnames(y)$x2[2:10],
        x3=dimnames(y)$x3[2:10],
        t =dimnames(y)$t,
        i =dimnames(y)$i,
        e =c("e1", "e2", "e3")
    )
)
for (x1 in 1:dim(cy)[1])
    for (x2 in 1:dim(cy)[2])
        for (x3 in 1:dim(cy)[3])
            for(t in 1:dim(cy)[4])
                for (i in 1:dim(cy)[5])
                    cy[x1, x2, x3, t, i, ] <- sort(eigen(ddy[x1, x2, x3, t, i, , ])$values)
cy %>% dim

cy[, 5, 5, 21, 3, ]

cy[5, , 5, 21, 3, ]

cy[5, 5, , 21, 3, ]

yt <- tucker(as.tensor(y), c(3, 3, 3, 3, 3))

yt

y %>% dim

midpoints <- function(x) (x[-1] + x[-length(x)]) / 2

y.scaled <- apply(y, 1:4, scale)
y.scaled %>% dim

y.deltas <- array(
    0,
    dim = c(
        dim(y.scaled)[1],
        dim(y.scaled)[2] - 1,
        dim(y.scaled)[3] - 1,
        dim(y.scaled)[4] - 1,
        dim(y.scaled)[5],
        3
    ),
    dimnames = list(
        i = dimnames(y)$i,
        x1 = midpoints(as.numeric(dimnames(y.scaled)$x1)),
        x2 = midpoints(as.numeric(dimnames(y.scaled)$x2)),
        x3 = midpoints(as.numeric(dimnames(y.scaled)$x3)),
        t = dimnames(y)$t,
        d = c("+++", "++-", "+--")
    )
)
im <- 1 : dim(y.deltas)[2]
ip <- 1 + im
y.deltas[, , , , , 1] <- (y.scaled[, ip, ip, ip, ] + y.scaled[, im, im, im, ]) / 2
y.deltas[, , , , , 2] <- (y.scaled[, ip, ip, im, ] + y.scaled[, im, im, ip, ]) / 2
y.deltas[, , , , , 3] <- (y.scaled[, ip, im, im, ] + y.scaled[, im, ip, ip, ]) / 2
y.deltas %>% dim

y.maxdelta <- apply(
    apply(y.deltas, 1:5, sd),
    2:4,
    max
)
y.maxdelta %>% dim

z <- data.table(melt(y.maxdelta))[, rank := rank(value, ties="random")]
z %>% head

ggplot(z, aes(x=x2, y=x3, fill=value)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_distiller(palette = "Spectral") +
    facet_grid(x1 ~ .)

ggplot(z, aes(x=x2, y=x3, fill=rank)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_distiller(palette = "Spectral") +
    facet_grid(x1 ~ .)

ggpairs(s3d, columns=c(2:4, 6:8), mapping=aes(color=factor(t)))
