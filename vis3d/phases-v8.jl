using Pkg
Pkg.activate(".")

using InfoVis
using InfoVis.Primitives

import CSV
import Colors: @colorant_str
import DataFrames: DataFrame
import InfoVis.Plots: axes, scatterplot
import InfoVis.Transport: connect, responses, stop
import StaticArrays: SVector
import StatsBase: UnitRangeTransform, fit, transform


function normalize(xs)
  dt = fit(UnitRangeTransform, xs, dims=1)
  transform(dt, xs)
end


z = connect(string("wss://", false ? "substrate.functionally.dev" : "127.0.0.1", ":42041/infovis/v4/julia"));

sleep(5);

z |> resetdisplay();

z |> axes(labels = ["x1", "x2", "x3"])

xs = CSV.read("../xs-3d-20200322a.csv");
ys = CSV.read("../ys-3d-20200322a.csv");

xys = join(xs, ys, on=:case);

v = DataFrame(
  id     = xys.case .* 1000 .+ round.(Int64, xys.t .* 10),
  x1     = xys.x1                                        ,
  x2     = xys.x2                                        ,
  x3     = xys.x3                                        ,
  f      = round.(Int64, 1 .+ 2 .* xys.t)                ,
  y1     = normalize(xys.y1)                             ,
  y2     = normalize(xys.y2)                             ,
  y3     = normalize(xys.y3)                             ,
);

α = 0.075

z |> request([
  Axis(
    row.id                                                                ,
    SVector(row.x1             , row.x2             , row.x3             ),
    SVector(row.x1 + α * row.y1, row.x2 + α * row.y2, row.x3 + α * row.y3),
    size  = 0.005                                                         ,
    color = colorant"yellow"                                              ,
  )
  for row in eachrow(v[v.f .== 21, :])
])

z |> rectangles(1, [SVector(
    SVector(0.587943644718903, 0., 0.),
    SVector(0.587943644718903, 0., 1.),
    SVector(0.587943644718903, 1., 0.),
)], color = colorant"#1b9e77ff")

z |> rectangles(2, [SVector(
    SVector(0., 0.955400602242447, 0.),
    SVector(0., 0.955400602242447, 1.),
    SVector(1., 0.955400602242447, 0.),
)], color = colorant"#d95f02ff")

z |> rectangles(3, [SVector(
    SVector(0.587943644718903, 0.451441434343024, 0.),
    SVector(0.587943644718903, 0.451441434343024, 1.),
    SVector(1., 0.451441434343024, 0.),
)], color = colorant"#7570b3ff")

z |> rectangles(4, [SVector(
    SVector(0., 0., 0.837110983745585),
    SVector(0., 1., 0.837110983745585),
    SVector(1., 0., 0.837110983745585),
)], color = colorant"#e7298aff")

z |> rectangles(5, [SVector(
    SVector(0., 0., 0.570905262263832),
    SVector(0., 1., 0.570905262263832),
    SVector(0.587943644718903, 0., 0.570905262263832),
)], color = colorant"#66a61eff")

z |> rectangles(6, [SVector(
    SVector(0.587943644718903, 0.451441434343024, 0.33455438805106),
    SVector(0.587943644718903, 0.955400602242447, 0.33455438805106),
    SVector(1., 0.451441434343024, 0.33455438805106),
)], color = colorant"darkblue")
