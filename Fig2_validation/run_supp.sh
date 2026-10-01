#!/bin/bash
# Regenerates supplementary/Supp_{CD177,Parasite}_pipelines.{pdf,svg,png}: the step-by-step pipeline tables
# that Fig. 2 replaced (steps only; times / verdicts are in the Fig. 2 scorecard). Needs no run data.
# Archive the previous outputs into archive/supp_vN_YYYY-MM-DD/ before re-running.
HERE="$(cd "$(dirname "$0")" && pwd)"
docker run --rm --user $(id -u):$(id -g) -e MPLCONFIGDIR=/tmp -v "$HERE":/work \
  --entrypoint /opt/conda/envs/local_imagent_J/bin/python agenticj:gpu-local /work/make_supp_pipelines.py
