
require(abind)
require(data.table)
require(magrittr)

require(ggplot2)

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

dy1 <- y[2:11, , , , ] - y[1:10, , , , ]
dy2 <- y[, 2:11, , , ] - y[, 1:10, , , ]
dy3 <- y[, , 2:11, , ] - y[, , 1:10, , ]

dy <- abind(
    (dy1[, 2:11, 2:11, , ] + dy1[, 1:10, 1:10, , ]) / 2,
    (dy2[2:11, , 2:11, , ] + dy2[1:10, , 1:10, , ]) / 2,
    (dy3[2:11, 2:11, , , ] + dy3[1:10, 1:10, , , ]) / 2,
    along=6
)
dy %>% dim

ddy1 <- dy[2:10, , , , , ] - dy[1:9, , , , , ]
ddy2 <- dy[, 2:10, , , , ] - dy[, 1:9, , , , ]
ddy3 <- dy[, , 2:10, , , ] - dy[, , 1:9, , , ]

ddy <- abind(
    (ddy1[, 2:10, 2:10, , , ] + ddy1[, 1:9, 1:9, , , ]) / 2,
    (ddy2[2:10, , 2:10, , , ] + ddy2[1:9, , 1:9, , , ]) / 2,
    (ddy3[2:10, 2:10, , , , ] + ddy3[1:9, 1:9, , , , ]) / 2,
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
        dim(ddy)[5]
    ),
    dimnames = list(
        x1=dimnames(y)$x1[2:10],
        x2=dimnames(y)$x2[2:10],
        x3=dimnames(y)$x3[2:10],
        t =dimnames(y)$t,
        i =dimnames(y)$i
    )
)
for (x1 in 1:dim(cy)[1])
    for (x2 in 1:dim(cy)[2])
        for (x3 in 1:dim(cy)[3])
            for(t in 1:dim(cy)[4])
                for (i in 1:dim(cy)[5])
                    cy[x1, x2, x3, t, i] <- det(ddy[x1, x2, x3, t, i, , ])
cy %>% dim

curv <- merge(
    s3d,
    dcast(
        data.table(melt(cy)),
        x1 + x2 + x3 + t ~ i
    )[, .(x1, x2, x3, t, c1=`1`, c2=`2`, c3=`3`)],
    by=c("x1", "x2", "x3", "t")
)[, .(case, t, x1, x2, x3, y1, y2, y3, c1, c2, c3)]
curv %>% head

save(xs3d, ys3d, s3d, y, dy, ddy, curv, file="visualizations-v4.rdata")

ggplot(curv, aes(x=x1, y=x2, color=t, size=c1)) + geom_point() + ggtitle("c1")

ggplot(curv, aes(x=x2, y=x3, color=t, size=c1)) + geom_point() + ggtitle("c1")

ggplot(curv, aes(x=x3, y=x1, color=t, size=c1)) + geom_point() + ggtitle("c1")

ggplot(curv, aes(x=x1, y=x2, color=t, size=c2)) + geom_point() + ggtitle("c2")

ggplot(curv, aes(x=x2, y=x3, color=t, size=c2)) + geom_point() + ggtitle("c2")

ggplot(curv, aes(x=x3, y=x1, color=t, size=c2)) + geom_point() + ggtitle("c2")

ggplot(curv, aes(x=x1, y=x2, color=t, size=c3)) + geom_point() + ggtitle("c3")

ggplot(curv, aes(x=x2, y=x3, color=t, size=c3)) + geom_point() + ggtitle("c3")

ggplot(curv, aes(x=x3, y=x1, color=t, size=c3)) + geom_point() + ggtitle("c3")

cy <- array(
    0,
    dim = c(
        dim(ddy)[1],
        dim(ddy)[2],
        dim(ddy)[3],
        dim(ddy)[4],
        dim(ddy)[5]
    ),
    dimnames = list(
        x1=dimnames(y)$x1[2:10],
        x2=dimnames(y)$x2[2:10],
        x3=dimnames(y)$x3[2:10],
        t =dimnames(y)$t,
        i =dimnames(y)$i
    )
)
for (x1 in 1:dim(cy)[1])
    for (x2 in 1:dim(cy)[2])
        for (x3 in 1:dim(cy)[3])
            for(t in 1:dim(cy)[4])
                for (i in 1:dim(cy)[5])
                    cy[x1, x2, x3, t, i] <- max(abs(ddy[x1, x2, x3, t, i, , ]))
cy %>% dim

curv <- merge(
    s3d,
    dcast(
        data.table(melt(cy)),
        x1 + x2 + x3 + t ~ i
    )[, .(x1, x2, x3, t, c1=`1`, c2=`2`, c3=`3`)],
    by=c("x1", "x2", "x3", "t")
)[, .(case, t, x1, x2, x3, y1, y2, y3, c1, c2, c3)]
curv %>% head

ggplot(curv, aes(x=x1, y=x2, color=t, size=c1)) + geom_point() + ggtitle("c1")

ggplot(curv, aes(x=x2, y=x3, color=t, size=c1)) + geom_point() + ggtitle("c1")

ggplot(curv, aes(x=x3, y=x1, color=t, size=c1)) + geom_point() + ggtitle("c1")

ggplot(curv, aes(x=x1, y=x2, color=t, size=c2)) + geom_point() + ggtitle("c2")

ggplot(curv, aes(x=x2, y=x3, color=t, size=c2)) + geom_point() + ggtitle("c2")

ggplot(curv, aes(x=x3, y=x1, color=t, size=c2)) + geom_point() + ggtitle("c2")

ggplot(curv, aes(x=x1, y=x2, color=t, size=c3)) + geom_point() + ggtitle("c3")

ggplot(curv, aes(x=x2, y=x3, color=t, size=c3)) + geom_point() + ggtitle("c3")

ggplot(curv, aes(x=x3, y=x1, color=t, size=c3)) + geom_point() + ggtitle("c3")
