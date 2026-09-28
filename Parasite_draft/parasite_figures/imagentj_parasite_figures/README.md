# Parasite-task figures — three figures, the scripts, and the data behind them

Three figures comparing how a published manual analysis and three
agent configurations performed the same three bioimage-analysis tasks on a
*Leishmania major* / macrophage dataset. Captions are in `figures/CAPTIONS.md`.

    figures/       the three PNGs, the same three as SVG, and the captions
    scripts/       everything needed to regenerate them
    data/          the inputs the scripts read — nothing else is required

Every number in every figure comes from a run's own output or from the article's
source data. Nothing is re-measured here.

## Regenerating the figures

Needs Python 3 with `numpy`, `scipy`, `matplotlib`, `scikit-image` and
`tifffile`. It does **not** need `openpyxl`, Java, ImageJ or Bio-Formats — the
workbooks and the ImageJ ROI files are read directly by `xlsx_min.py` and
`ij_roi.py`.

    cd scripts

    python3 figure_task1_overlays.py \
        --gt-root ../data/parasite_task_results \
        --one-frame ../data/one_frame \
        --out ../figures/figure_task1_overlays.png

    python3 figure_tasks23_dotplots.py \
        --gt-root ../data/parasite_task_results \
        --published ../data/published/jci.insight.169020.sdval.xlsx \
        --out ../figures/figure_tasks23_dotplots.png

    python3 figure_workflows.py --out-dir ../figures --tasks 1

Each prints the values it plotted, so a rebuild can be checked against the
figure rather than trusted. `figure_workflows.py` reads no data at all: its
content is the runs' recorded pipelines, written into the script.

## What the reference is

The laboratory's analysis, in the ImageJ ROI sets and Excel workbooks it was
done in. The proliferation index throughout is

    corrected green = parasite green - mean(3 local background regions, green)
    corrected red   = parasite red   - mean(3 local background regions, red)
    index           = 10 - corrected red / corrected green

`data/published/jci.insight.169020.sdval.xlsx` is the source data of JCI Insight
169020, whose sheets `4C` and `4D` are the published figures. Those values agree
with an independent parse of the laboratory's workbooks to the last digit — 6
values for task 1, 15 for task 2, all 828 for task 3 — which is why one column
can be labelled "published" and used as the reference at the same time.

Two properties of those workbooks are worth knowing before reading any of the
numbers:

* `Control<n>` is a parasite in a TUNEL-**negative** host. It is not a background
  region; background regions appear as `Control area <n>` elsewhere in the same
  sheet.
* Task 2's four workbooks describe **15** parasites, not 36. Each holds two
  apoptotic parasites unique to its subset plus the same seven controls, copied
  into all four.

## Contents of data/

    parasite_task_results/
      oneframe_results/          the laboratory's task-1 workbook
      auto-mode/task_1,2/        the Advanced configuration's outputs
      fast-mode/.../task_1,2/    the Fast configuration's outputs
      hil_results/repeat_1,2/    the two agent-plus-human repeats
    one_frame/                   the ImageJ ROI files: the six nuclei, their
                                 backgrounds, and the overview rectangle
    published/                   the article's source data

The four TIFF images were re-saved with lossless deflate compression, 34 MB down
to 8 MB, and verified pixel-identical to the originals. Only the files these
three figures read are included; the full study is larger.
