import itertools as it
import numpy     as np
import pandas    as pd

from math import nan


class Sampler:

    f  = None
    xs = None
    ts = None
    ys = None

    _digits   = 30
    _scale    = 1 << _digits
    _maxdepth = - np.int(np.log2(np.sqrt(np.finfo(float).eps)))

    def _level(self, i):
        if (i == 0):
            return self._digits
        n = 0
        while i & 1 == 0:
            n = n + 1
            i = i >> 1
        return n

    def _depth(self, i):
        return self._digits - self._level(i)

    def _corner(self, i, offset = 0):
        d = self._level(i) - offset
        return [i - (1 << d), i + (1 << d)]

    def _corners(self, i1, i2, i3, offset = 0):
        def axis(s1, s2, s3):
            if s1 == 0:
                return 4 * s1 + 2 * s2 + s3
            else:
                return 4 * (1 - s1) + 2 * (1 - s2) + (1 - s3)
        result = pd.DataFrame(
            {"d1": [1], "d2" : [2], "d3" : [4]}
        ).merge(
            pd.DataFrame({"d1" : 1, "i1" : self._corner(i1, offset=offset), "s1" : [0, 1]})
        ).merge(
            pd.DataFrame({"d2" : 2, "i2" : self._corner(i2, offset=offset), "s2" : [0, 1]})
        ).merge(
            pd.DataFrame({"d3" : 4, "i3" : self._corner(i3, offset=offset), "s3" : [0, 1]})
        )
        result["axis"] = result.apply(lambda row: axis(row.s1, row.s2, row.s3), axis=1)
        return result[["i1", "i2", "i3", "axis"]].set_index(["i1", "i2", "i3"])

    def _candidates(self, i1, i2, i3):
        result = pd.DataFrame(np.array(list(it.product(*[
            [i1] + self._corner(i1),
            [i2] + self._corner(i2),
            [i3] + self._corner(i3),
        ]))), columns=["i1", "i2", "i3"])
        result["x1"] = result["i1"] / self._scale
        result["x2"] = result["i2"] / self._scale
        result["x3"] = result["i3"] / self._scale
        return result.set_index(["i1", "i2", "i3"])

    def __init__(self, f, ts):
        self.f = f
        self.ts = ts
        self.xs = pd.DataFrame(np.array(list(it.product(*[
            [0, self._scale],
            [0, self._scale],
            [0, self._scale],
        ]))), columns=["i1", "i2", "i3"])
        self.xs = self.xs.assign(**{
            "sequence"  : range(1, 9)                ,
            "generation": 0                          ,
            "x1"        : self.xs["i1"] / self._scale,
            "x2"        : self.xs["i2"] / self._scale,
            "x3"        : self.xs["i3"] / self._scale,
            "compute"   : True                       ,
            "measure"   : False                      ,
            "probed"    : True                       ,
            "s1"        : 0                          ,
            "s2"        : 0                          ,
            "s3"        : 0                          ,
            "s"         : 0                          ,
        }).append(
            pd.DataFrame({
                "i1"         : [self._scale >> 1],
                "i2"         : [self._scale >> 1],
                "i3"         : [self._scale >> 1],
                "sequence"   : [9]               ,
                "generation" : [-1]              ,
                "x1"         : [0.5]             ,
                "x2"         : [0.5]             ,
                "x3"         : [0.5]             ,
                "compute"    : [False]           ,
                "measure"    : [True]            ,
                "probed"     : [False]           ,
                "s1"         : [nan]             ,
                "s2"         : [nan]             ,
                "s3"         : [nan]             ,
                "s"          : [nan]             ,
            })
        ).set_index(["i1", "i2", "i3"])
        self.ys = None

    def compute(self):
        needed = self.xs[self.xs.compute]
        for row in needed.itertuples():
            try:
                fx = self.f([row.x1, row.x2, row.x3], self.ts)
            except Exception:
                print("Exception evaluating function at " + str(row) + ".")
                fx = []
            y = pd.DataFrame(
                fx if len(fx) == len(self.ts) else np.full((len(self.ts), 3), nan),
                columns=["y1", "y2", "y3"]
            )
            y["sequence"] = row.sequence
            y["t"       ] = self.ts
            y = y.set_index(["sequence"])[["t", "y1", "y2", "y3"]]
            self.ys = y.append(self.ys)
        self.xs.sort_index(inplace=True)
        self.xs.loc[needed.index, "compute"] = False

    def measure(self, focus = 2):
        result = None
        for row in self.xs[self.xs.measure].itertuples():
            corners = self._corners(row.Index[0], row.Index[1], row.Index[2])
            corners["center"] = row.sequence
            result = corners.append(result)
        result = result.join(
            self.xs
        ).set_index(
            ["sequence"]
        )[
            ["center", "axis"]
        ].join(
            self.ys
        ).groupby(
            ["center", "axis", "t"]
        ).agg(
            np.mean
        ).reset_index(
        ).drop(
            columns=["t"]
        ).groupby(
            ["center", "axis"]
        ).agg(
            np.std
        ).reset_index(
        ).drop(
            columns=["axis"]
        ).groupby(
            ["center"]
        ).agg(
            np.max
        ).rename_axis(
            index={"center" : "sequence"}
        ).rename(
            columns={"y1" : "s1", "y2" : "s2", "y3" : "s3"}
        )

        result["generation"] = np.max(self.xs.generation) + 1

        self.xs.reset_index(inplace=True)
        self.xs.set_index(["sequence"], inplace=True)
        self.xs.update(result)
        self.xs.reset_index(inplace=True)

        minmax = self.ys.drop(columns=["t"]).aggregate([np.min, np.max])
        def normalize(i, s1, s2, s3):
            return max(
                s1 / (minmax.iloc[1, 0] - minmax.iloc[0, 0]),
                s2 / (minmax.iloc[1, 1] - minmax.iloc[0, 1]),
                s3 / (minmax.iloc[1, 2] - minmax.iloc[0, 2]),
            ) * focus**self._depth(i)

        self.xs["measure"] = False
        self.xs["s"] = self.xs.apply(
            lambda row: normalize(row["i1"], row["s1"], row["s2"], row["s3"]),
            axis=1
        )
        self.xs.set_index(["i1", "i2", "i3"], inplace=True)
        self.xs.sort_index(inplace=True)

    def probe(self, alpha = 1, maxdepth = _maxdepth):
        choice = self.xs[
            ~self.xs.probed & ~np.isnan(self.xs.s)
        ][
            ["sequence", "s"]
        ].sample(1, weights="s")
        candidates = self._candidates(choice.iloc[0].name[0], choice.iloc[0].name[1], choice.iloc[0].name[2])
        candidates = candidates.loc[candidates.index.difference(self.xs.index)]
        n = max(self.xs.sequence)
        if len(candidates) > 0:
            candidates = candidates.assign(**{
                "sequence"   : range(n + 1, n + len(candidates) + 1),
                "generation" : -1                                   ,
                "compute"    : True                                 ,
                "measure"    : False                                ,
                "probed"     : False                                ,
                "s1"         : nan                                  ,
                "s2"         : nan                                  ,
                "s3"         : nan                                  ,
                "s"          : nan                                  ,
            })
            n = max(candidates.sequence)
        probes = self._corners(choice.iloc[0].name[0], choice.iloc[0].name[1], choice.iloc[0].name[2], offset=1)
        probes = probes.assign(**{
            "x1"         : probes.index.get_level_values("i1").values / self._scale,
            "x2"         : probes.index.get_level_values("i2").values / self._scale,
            "x3"         : probes.index.get_level_values("i3").values / self._scale,
            "sequence"   : range(n + 1, n + len(probes) + 1)                       ,
            "generation" : -1                                                      ,
            "compute"    : False                                                   ,
            "measure"    : True                                                    ,
            "probed"     : False                                                   ,
            "s1"         : nan                                                     ,
            "s2"         : nan                                                     ,
            "s3"         : nan                                                     ,
            "s"          : nan                                                     ,
        })
        result = self.xs.append(probes[["sequence", "generation", "x1", "x2", "x3", "compute", "measure", "probed", "s1", "s2", "s3", "s"]]).sort_index()
        if len(candidates) > 0:
            result = result.append(candidates[["sequence", "generation", "x1", "x2", "x3", "compute", "measure", "probed", "s1", "s2", "s3", "s"]]).sort_index()
        result.loc[result.sequence == choice.sequence[0], ["compute", "measure", "probed"]] = [True, False, True]
        self.xs = result

    def cycle(self, n = 1, focus = 2, alpha = 1, maxdepth = _maxdepth):
        for i in range(n):
            self.compute()
            self.measure(focus = focus)
            self.probe(alpha = alpha, maxdepth = maxdepth)
        return (self.xs.shape, self.ys.shape)
