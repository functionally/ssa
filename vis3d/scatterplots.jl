using InfoVis
using InfoVis.Primitives

import CSV
import Colors: colormap
import DataFrames: DataFrame
import InfoVis.Plots: axes, scatterplot
import InfoVis.Transport: connect, responses, stop
import StaticArrays: SVector
import StatsBase: UnitRangeTransform, fit, transform


function normalize(xs)
  dt = fit(UnitRangeTransform, xs, dims=1)
  transform(dt, xs)
end

cm = colormap("RdBu") |> reverse

function makecolor(x)
  cm[ceil(Int64, 1 + (length(cm) - 1) * x)]
end


z = connect("ws://104.198.152.159:42042/infovis/v4/julia");

@async responses(z)


xs = CSV.read("../xs-3d-20200322a.csv");
ys = CSV.read("../ys-3d-20200322a.csv");

xys = join(xs, ys, on=:case);

w = DataFrame(
  case   = xys.case                                      ,
  x1     = xys.x1                                        ,
  x2     = xys.x2                                        ,
  x3     = xys.x3                                        ,
  t      = xys.t                                         ,
  y1     = normalize(xys.y1)                             ,
  y2     = normalize(xys.y2)                             ,
  y3     = normalize(xys.y3)                             ,
  tcolor = makecolor.(xys.t ./ 10)                       ,
  id     = xys.case .* 1000 .+ round.(Int64, xys.t .* 10),
)


z |> resetdisplay();

f = 0;

cases = [
  (:y1, :y3, :y2),
  (:x1, :y2, :y1),
  (:x1, :y3, :y2),
  (:x1, :y1, :y3),
  (:x2, :y2, :y1),
  (:x2, :y3, :y2),
  (:x2, :y1, :y3),
  (:x3, :y2, :y1),
  (:x3, :y3, :y2),
  (:x3, :y1, :y3),
  (:x1, :y1, :x2),
  (:x1, :y2, :x2),
  (:x1, :y3, :x2),
  (:x2, :y1, :x3),
  (:x2, :y2, :x3),
  (:x2, :y3, :x3),
  (:x3, :y1, :x1),
  (:x3, :y2, :x1),
  (:x3, :y3, :x1),
];

for (a, b, c) in cases
  global f += 1
  z |> axes(labels = [string(a), string(b), string(c)], frame = f) |>
       scatterplot(w, :id, a, b, c, color = :tcolor, frame = f)
end

z |> stop
