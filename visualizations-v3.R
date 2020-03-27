
require(data.table)
require(magrittr)

require(ggplot2)
require(GGally)

xs2d <- fread("xs-2d-20200322a.csv")[, case:=factor(case)]
xs2d %>% dim

ys2d <- fread("ys-2d-20200322a.csv")[, case:=factor(case)]
ys2d %>% dim

s2d = merge(xs2d, ys2d)
s2d %>% head

xs3d <- fread("xs-3d-20200322a.csv")[, case:=factor(case)]
xs3d %>% dim

ys3d <- fread("ys-3d-20200322a.csv")[, case:=factor(case)]
ys3d %>% dim

s3d = merge(xs3d, ys3d)
s3d %>% head

jitter2d <- function(df)
    df[, .(
        case,
        x1=jitter(x1, factor=1),
        x2=jitter(x2, factor=1),
        t,
        y1=jitter(y1, factor=1),
        y2=jitter(y2, factor=1)
    )]

jitter3d <- function(df)
    df[, .(
        case,
        x1=jitter(x1, factor=1),
        x2=jitter(x2, factor=1),
        x3=jitter(x3, factor=1),
        t,
        y1=jitter(y1, factor=1),
        y2=jitter(y2, factor=1),
        y3=jitter(y3, factor=1)
    )]

levels2 <- c(0, 1)
levels3 <- c(0, 0.5, 1)
levels4 <- c(0, 0.3, 0.7, 1)
levels5 <- c(0, 0.2, 0.5, 0.8, 1)

ggparcoord(
    s2d[x1 %in% levels2 & x2 %in% levels3 & t == 10] %>% jitter2d,
    columns=c(2:3,5:6),
    groupColumn=1
)

ggparcoord(
    s2d[t == 10] %>% jitter2d,
    columns=c(2:3,5:6),
    groupColumn=1
) + guides(color=FALSE)

ggparcoord(
    s2d[t == 5] %>% jitter2d,
    columns=c(2:3,5:6),
    groupColumn=1
) + guides(color=FALSE)

ggparcoord(
    s3d[x1 %in% levels3 & x2 %in% levels3 & x3 %in% levels3 & t == 10] %>% jitter3d,
    columns=c(2:4,6:8),
    groupColumn=1
)

ggparcoord(
    s3d[t == 10] %>% jitter3d,
    columns=c(2:4,6:8),
    groupColumn=1
) + guides(color=FALSE)

ggparcoord(
    s3d[t == 5] %>% jitter3d,
    columns=c(2:4,6:8),
    groupColumn=1
) + guides(color=FALSE)

ggpairs(
    s3d[t == 10],
    columns=c(2:4,6:8)
)

ggpairs(
    s3d[t == 10],
    columns=6:8
)

ggpairs(s3d, columns=6:8, mapping=aes(color=factor(t)))

x2.lab <- paste("x2 =", levels4)
names(x2.lab) <- levels4
x3.lab <- paste("x3 =", levels4)
names(x3.lab) <- levels4

for (x1.value in levels4) {
    g <- ggplot(
        melt(
            s3d[x1 == x1.value & x2 %in% levels4 & x3 %in% levels4],
            id.vars=c("case", "t", "x1", "x2", "x3")
        ),
        aes(x=t, y=value, color=variable)
    ) +
        facet_grid(
            x3 ~ x2,
            labeller=labeller(x2=x2.lab, x3=x3.lab)) +
        ggtitle(paste("x1 =", x1.value)) +
        geom_line()
    print(g)
}
