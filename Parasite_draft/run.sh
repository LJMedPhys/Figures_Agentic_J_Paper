#!/bin/bash
# Regenerates Fig_Parasite.{pdf,svg,png} from the bundled run outputs + article source data (parasite_figures/).
HERE="$(cd "$(dirname "$0")" && pwd)"
docker run --rm --user $(id -u):$(id -g) -e MPLCONFIGDIR=/tmp \
  -v "$HERE/parasite_figures/imagentj_parasite_figures":/p:ro -v "$HERE":/work \
  --entrypoint /opt/conda/envs/local_imagent_J/bin/python agenticj:gpu-local /work/make_figure.py
