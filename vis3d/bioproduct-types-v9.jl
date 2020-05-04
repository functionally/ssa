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


z = connect(string("wss://", true ? "substrate.functionally.dev" : "127.0.0.1", ":42041/infovis/v4/julia"));

sleep(5);


xs = CSV.read("../framework/xs-v9.tsv");

v = DataFrame(
  id     = xs.sequence,
  x1     = xs.x1      ,
  x2     = xs.x2      ,
  x3     = xs.x3      ,
);

z |> resetdisplay();

z |> axes(labels = ["bioproduct actual cost", "consumer inertia", "retrofit delay"]) |>
     scatterplot(v, :id, :x1, :x2, :x3, size = 0.0025, color = colorant"seagreen")

z |> sleep(60)

z |> stop
