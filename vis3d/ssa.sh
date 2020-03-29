#!/usr/bin/env bash

set -x

../../go-infovis --stderrthreshold=WARNING                    \
                 --demo 0.0.0.0:42041 /infovis/v4/julia-demo  \
                 ssa-*.pbb
