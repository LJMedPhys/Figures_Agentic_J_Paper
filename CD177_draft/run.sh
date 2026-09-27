#!/bin/bash
# Regenerates Fig_CD177.{pdf,svg,png} from the run outputs (read-only mounts) and inputs/.
# Helper scripts: agreement.py -> inputs/mask_agreement_per_cell.csv (agent vs published masks),
#                 published_from_crops.py -> inputs/published_cells_from_crops.csv (paper's own crops, s3 rules).
HERE="$(cd "$(dirname "$0")" && pwd)"
Y=/mnt/eternus/users/Yu/project/bioimage_agent; SUB=$Y/bioimage-agent-bench-run/outputs/submissions
AR=$SUB/agentic_j_cloud_gpt56/run_20260921_105141/cd177-neutrophil-morphology-raw-only/artifacts
FR=$SUB/agentic_j_cloud_gpt56_quick/run_20260924_130802/cd177-neutrophil-morphology-raw-only/artifacts
CR=$SUB/codex_cli_cloud_gpt56/run_20260921_192946/cd177-neutrophil-morphology-raw-only/artifacts
HP=$Y/Imagent_J/data/projects/cd177_neutrophil_morphology_microsam/data
docker run --rm --user $(id -u):$(id -g) -e MPLCONFIGDIR=/tmp \
  -v $AR/cell_features:/auto_cells:ro -v $FR/feature_statistics:/fast_stats:ro -v $FR/cell_features:/fast_cells:ro \
  -v $CR/feature_statistics:/codex_stats:ro -v $CR/cell_features:/codex_cells:ro \
  -v $HP/statistics:/hil_stats:ro -v $HP/measurements:/hil_cells:ro -v "$HERE":/work \
  --entrypoint /opt/conda/envs/local_imagent_J/bin/python agenticj:gpu-local /work/make_figure.py
