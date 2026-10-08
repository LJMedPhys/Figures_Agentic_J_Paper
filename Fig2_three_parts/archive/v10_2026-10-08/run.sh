#!/bin/bash
# Regenerates Fig2_three_parts_{orientationj,coloc2}.{pdf,svg,png}: Fig. 2 as (1) plugin cases, (2) CD177, (3) Leishmania.
# Reads the draft, Fig2_validation, BoneJ and Coloc 2 folders read-only (they are not modified). Archive the previous outputs into
# archive/vN_YYYY-MM-DD/ before re-running.
HERE="$(cd "$(dirname "$0")" && pwd)"
FIGS="$(dirname "$HERE")"
# CD177 data (all runs, + human = HIL6) is read from inputs/cd177_zenodo (Zenodo record), no external mounts
# parasite, Conventional Agentic-J + human: repeat 3
HIL3=/home/marilin/Documents/imagej_automation/parasite_task_results/hil_results/repeat_3/leishmania_tunel_position40
IMG=agenticj:gpu-local; docker image inspect $IMG >/dev/null 2>&1 || IMG=agenticj:cpu-local-napari
for V in orientationj coloc2; do
docker run --rm --user $(id -u):$(id -g) -e MPLCONFIGDIR=/tmp \
  -v "$FIGS/CD177_draft":/cd177:ro \
  -v "$FIGS/Parasite_draft/parasite_figures/imagentj_parasite_figures":/p:ro \
  -v "$FIGS/Parasite_draft/inputs":/parin:ro \
  -v "$HIL3":/hil3:ro \
  -v "$FIGS/Fig2_validation/inputs":/fv:ro -v "$FIGS/BoneJ_reproduction/inputs":/bonej:ro -v "$FIGS/Coloc2_reproduction/inputs":/coloc:ro \
  -v "$HERE":/work -e FIG2_VARIANT="$V" \
  --entrypoint /opt/conda/envs/local_imagent_J/bin/python $IMG /work/make_fig2_three.py
done
