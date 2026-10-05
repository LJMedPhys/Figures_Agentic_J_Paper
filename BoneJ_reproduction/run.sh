#!/bin/bash
# Regenerates Fig_BoneJ_reproduction.{pdf,svg,png} from inputs/ (copies of the BoneJ scoring files).
# Archive the previous outputs into archive/vN_YYYY-MM-DD/ before re-running.
HERE="$(cd "$(dirname "$0")" && pwd)"
docker run --rm --user $(id -u):$(id -g) -e MPLCONFIGDIR=/tmp -v "$HERE":/work \
  --entrypoint /opt/conda/envs/local_imagent_J/bin/python agenticj:gpu-local /work/make_fig_bonej.py
