#!/usr/bin/env bash

set -x

../../go-infovis --stderrthreshold=WARNING               \
                 --demo 0.0.0.0:42043 /infovis/v4/client \
                 ../certificate.crt                      \
                 ../private.key                          \
                 ssa-v5-*.pbb
