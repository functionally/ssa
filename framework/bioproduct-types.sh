#!/usr/bin/env nix-shell
#!nix-shell -i bash -p pythonEnv wineWowPackages.minimal

Xvfb :9 -screen 0 1024x768x16 &
XVFB_PID=$!

export DISPLAY=:9.0

python bioproduct-types.py

kill $XVFB_PID
