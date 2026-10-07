#!/bin/bash
# Regenerates Fig2_three_parts_{orientationj,coloc2}.{pdf,svg,png}: Fig. 2 as (1) plugin cases, (2) CD177, (3) Leishmania.
# Reads the draft, Fig2_validation, BoneJ and Coloc 2 folders read-only (they are not modified). Archive the previous outputs into
# archive/vN_YYYY-MM-DD/ before re-running.
HERE="$(cd "$(dirname "$0")" && pwd)"
FIGS="$(dirname "$HERE")"
# CD177 run outputs (panel d), same sources as ../CD177_draft/run.sh
Y=/mnt/eternus/users/Yu/project/bioimage_agent; SUB=$Y/bioimage-agent-bench-run/outputs/submissions
AR=$SUB/agentic_j_cloud_gpt56/run_20260921_105141/cd177-neutrophil-morphology-raw-only/artifacts
FR=$SUB/agentic_j_cloud_gpt56_quick/run_20260924_130802/cd177-neutrophil-morphology-raw-only/artifacts
CR=$SUB/codex_cli_cloud_gpt56/run_20260921_192946/cd177-neutrophil-morphology-raw-only/artifacts
HP=$Y/Imagent_J/data/CD177_task/runs/20260928_213951_hil4/app_data/projects/cd177_neutrophil_morphology/data
# parasite, Conventional Agentic-J + human: repeat 3
HIL3=/home/marilin/Documents/imagej_automation/parasite_task_results/hil_results/repeat_3/leishmania_tunel_position40
# CD177 run outputs only exist on the eternus machine; without them make_fig2_three.py reuses panels f-g from archive/v5
CDM=(); [ -d "$AR" ] && CDM=(-v $AR/cell_features:/auto_cells:ro -v $FR/cell_features:/fast_cells:ro -v $CR/cell_features:/codex_cells:ro -v $HP:/hil_cells:ro)
IMG=agenticj:gpu-local; docker image inspect $IMG >/dev/null 2>&1 || IMG=agenticj:cpu-local-napari
for V in orientationj coloc2; do
docker run --rm --user $(id -u):$(id -g) -e MPLCONFIGDIR=/tmp \
  -v "$FIGS/CD177_draft":/cd177:ro \
  -v "$FIGS/Parasite_draft/parasite_figures/imagentj_parasite_figures":/p:ro \
  -v "$FIGS/Parasite_draft/inputs":/parin:ro \
  "${CDM[@]}" -v "$HIL3":/hil3:ro \
  -v "$FIGS/Fig2_validation/inputs":/fv:ro -v "$FIGS/BoneJ_reproduction/inputs":/bonej:ro -v "$FIGS/Coloc2_reproduction/inputs":/coloc:ro \
  -v "$HERE":/work -e FIG2_VARIANT="$V" \
  --entrypoint /opt/conda/envs/local_imagent_J/bin/python $IMG /work/make_fig2_three.py
done
