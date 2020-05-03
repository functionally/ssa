import numpy   as np
import pandas  as pd
import Sampler as sa
import Vensim  as vs


v = vs.Vensim(
    "simulation",
     2.50,
     mdl_name="bioproduct-types.mdl",
     executable="wine vensim/vensimdp.exe"
)

v.make_lst([
    "bioproduct actual cost",
    "consumer inertia"      ,
    "retrofit delay"        ,
    "Adopters"              ,
    "NonAdopters"           ,
    "Potential Adopters"    ,
])

v.make_cmd()


xrange = np.array([[5, 35], [0.1, 5], [1, 35]])

def bioproducts(xs, ts):
    v.make_cin({
        "bioproduct actual cost" : xrange[0, 0] * (1 - xs[0]) + xrange[0, 1] * xs[0],
        "consumer inertia"       : xrange[1, 0] * (1 - xs[1]) + xrange[1, 1] * xs[1],
        "retrofit delay"         : xrange[2, 0] * (1 - xs[2]) + xrange[2, 1] * xs[2],
    })
    v.make_cmd()
    v.run_vensim()
    return v.read_tsv()[["Adopters", "NonAdopters", "Potential Adopters"]].to_numpy()


z = sa.Sampler(bioproducts, np.arange(2015, 2051, 2.5))

while True:
    z.cycle(5, focus = 0.5, alpha = 0.5)
    z.xs.assign(**{
        "bioproduct actual cost" : xrange[0, 0] * (1 - z.xs.x1) + xrange[0, 1] * z.xs.x1,
        "consumer inertia"       : xrange[1, 0] * (1 - z.xs.x2) + xrange[1, 1] * z.xs.x2,
        "retrofit delay"         : xrange[2, 0] * (1 - z.xs.x3) + xrange[2, 1] * z.xs.x3,
    }).to_csv("xs.tsv", sep = "\t")
    z.ys.rename(columns={
        "t"  : "Year"              ,
        "y1" : "Adopters"          ,
        "y2" : "NonAdopters"       ,
        "y3" : "Potential Adopters",
    }).to_csv("ys.tsv", sep = "\t")

v.delete_run(delete_lst=True)
