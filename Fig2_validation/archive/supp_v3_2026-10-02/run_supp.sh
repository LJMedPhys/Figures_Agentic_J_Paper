#!/bin/bash
# Regenerates supplementary/Supp_{CD177,Parasite}_pipelines.{pdf,svg,png}: the step-by-step pipeline tables
# that Fig. 2 replaced, plus (parasite, panel b) the count-comparison table. Panel b reads the parasite run
# data through make_fig2.py's data block, hence the two parasite mounts (read-only).
# Archive the previous outputs into archive/supp_vN_YYYY-MM-DD/ before re-running.
HERE="$(cd "$(dirname "$0")" && pwd)"
FIGS="$(dirname "$HERE")"
docker run --rm --user $(id -u):$(id -g) -e MPLCONFIGDIR=/tmp \
  -v "$FIGS/Parasite_draft/parasite_figures/imagentj_parasite_figures":/p:ro \
  -v "$FIGS/Parasite_draft/inputs":/parin:ro \
  -v "$HERE":/work \
  --entrypoint /opt/conda/envs/local_imagent_J/bin/python agenticj:gpu-local /work/make_supp_pipelines.py
