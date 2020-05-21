
require(abind)
require(data.table)
require(deSolve)
require(LaplacesDemon)
require(magrittr)
require(np)
require(rTensor)

require(ggplot2)
require(GGally)

xs = fread(cmd="xzcat xs-v9.tsv.xz")
xs %>% dim

xs[!is.na(s)] %>% dim

ys = fread(cmd="xzcat ys-v9.tsv.xz")
ys %>% summary

xys = merge(xs, ys, by="sequence")[, .(
    `sequence`,
    `bioproduct actual cost`,
    `consumer inertia`,
    `retrofit delay`,
    `Year`,
    `NonAdopters`,
    `Potential Adopters`,
    `Adopters`
)]

ggpairs(
    xs,
    16:18,
#   mapping = aes(color=!probed),
    lower = list(continuous = wrap("points", alpha = 0.075, size=0.1))
)

Modes(xs[!is.na(s), `bioproduct actual cost`])$modes

Modes(xs[!is.na(s), `consumer inertia`])$modes

Modes(xs[!is.na(s), `retrofit delay`])$modes

ggplot(xs[!is.na(s)], aes(x=s)) + geom_histogram(bins=30)

ggplot(
    melt(xys, id.vars=c("sequence", "Year", "bioproduct actual cost", "consumer inertia", "retrofit delay")),
    aes(x=Year, y=value, group=sequence, color=`bioproduct actual cost`)) +
geom_line(alpha=0.125, size=0.25) +
scale_color_distiller(type="div") +
facet_grid(variable ~ .)

ggplot(
    melt(xys, id.vars=c("sequence", "Year", "bioproduct actual cost", "consumer inertia", "retrofit delay")),
    aes(x=Year, y=value, group=sequence, color=`consumer inertia`)) +
geom_line(alpha=0.125, size=0.25) +
scale_color_distiller(type="div") +
facet_grid(variable ~ .)

ggplot(
    melt(xys, id.vars=c("sequence", "Year", "bioproduct actual cost", "consumer inertia", "retrofit delay")),
    aes(x=Year, y=value, group=sequence, color=`retrofit delay`)) +
geom_line(alpha=0.125, size=0.25) +
scale_color_distiller(type="div") +
facet_grid(variable ~ .)

ggpairs(
    xys[Year == 2050],
    6:8,
#   mapping = aes(color=!probed),
    lower = list(continuous = wrap("points", alpha = 0.1, size=0.1))
)
