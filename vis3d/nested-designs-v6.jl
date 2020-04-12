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


z = connect(string("wss://", true ? "substrate.functionally.dev" : "127.0.0.1", ":42041/infovis/v4/julia"));

sleep(5);


ds = CSV.read("../xs-v6.csv");

u = DataFrame(
  id   = ds.sequence,
  x1   = ds.x1      ,
  x2   = ds.x2      ,
  x3   = ds.x3      ,
  size = 0.01       ,
)


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

z |> resetdisplay();

for f in unique(v.f)[end:end]
  ff = 21 - f
  z |> label(
    10000000 * f,
    string("t = ", (f - 1.) / 2),
    SVector(0., 1. , 0.),
    SVector(1., 1. , 0.),
    SVector(0., 1.1, 0.),
    size  = 0.1         ,
    frame = ff          ,
  )
  z |> axes(labels = ["x1", "x2", "x3"], frame = ff) |>
       scatterplot(u, :id, :x1, :x2, :x3, size = :size, color = colorant"seagreen", frame = ff)
  z |> request([
    Axis(
      row.id                                                                ,
      SVector(row.x1             , row.x2             , row.x3             ),
      SVector(row.x1 + α * row.y1, row.x2 + α * row.y2, row.x3 + α * row.y3),
      size  = 0.005                                                         ,
      color = colorant"orange"                                              ,
      frame = ff                                                            ,
    )
    for row in eachrow(v[v.f .== f, :])
  ])
end

z |> sleep(60)

z |> stop
