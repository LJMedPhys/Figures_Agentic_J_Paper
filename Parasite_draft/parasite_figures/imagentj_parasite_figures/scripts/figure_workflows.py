#!/usr/bin/env python3
"""How each configuration reached its answer: the pipeline, step by step.

    figure_workflows.py --out-dir .        # writes figure_workflow_task{1,2,3}.svg/.png

One figure per task. Each column is a configuration, each block a step, and the
small grey line under a step names the tool or method that performed it. The
badge under the configuration name is time to result; the line at the foot is
what came out.

WHERE THE STEPS COME FROM
-------------------------
Not from this script's imagination. Every block is taken from the run's own
record: `state_ledger.json` for the plan and the decisions, the QC report beside
each output for the rule actually applied, and the saved scripts where a run kept
them (see trace_methods.py, which mines all three). The manual column is the
laboratory's own description of how the published figures were made.

TIME TO RESULT
--------------
Advanced ran each task as its own session, so its times are per task as recorded.
Fast ran all three in ONE session and each result.json carries the cumulative
session usage, so its per-task times are differences: 22, then 91-22 = 69, then
91-91 = 0 -- task 3 completed no query before the deadline. Both configurations
are shown against the same 85 minute budget, so neither is credited with the
longer session it happened to be given.

The agent+human column needs more care. Repeat 1's usage log was never copied
out -- the file in its folder belongs to a different session and references none
of its outputs -- so repeat 1 has no per-task time. Repeat 2's does.

For TASK 1 that costs nothing, because the two repeats ran the same pipeline on
channel 2 and produced the same six corrected values to thirteen digits; the
column is therefore labelled "example repeat" and carries repeat 2's 32 minutes.
(Repeat 2 additionally ran channel 3 first and published that instead, which is
a difference in what it adopted, not in the steps drawn here.)

For TASKS 2 AND 3 the repeats genuinely diverge -- how repeated controls were
collapsed, and whether backgrounds were paired per host or by clustering -- so
those columns name repeat 1 explicitly and leave the time blank rather than
borrow a number from a run that did something else.
"""

import argparse
import os
import re
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

INK = "#1b1b1b"
MODES = {
    "Published (manual)":  dict(edge="#5b6770", fill="#eceef0"),
    "Agent + human":       dict(edge="#3d6d8f", fill="#e8eff5"),
    "Advanced":            dict(edge="#2a9d8f", fill="#e5f4f1"),
    "Fast":                dict(edge="#e07a1f", fill="#fdf0e2"),
}

TASKS = {
    1: dict(
        title="Task 1 — TUNEL classification of six nuclei, position 40",
        columns=[
            ("Published (manual)", "not recorded", [
                ("Draw ROIs", "by hand: 6 nuclei + 8 background regions"),
                ("Measure", "ImageJ ROI Manager, TUNEL channel"),
                ("Export", "ROI-set intensities to Excel"),
                ("Correct", "raw - mean(background), in Excel"),
                ("Classify", "TUNEL+ by inspection"),
            ], "6 nuclei, 2 called TUNEL+"),
            ("Agent + human (example repeat)", "32 min", [
                ("Import", "Bio-Formats windowless, series 39, native 16-bit"),
                ("Place the given ROIs", "RoiDecoder + overview-origin translation"),
                ("Measure", "ImageJ ImageStatistics, no rescaling"),
                ("Correct", "unweighted mean of the 8 given regions"),
                ("Classify", "largest adjacent gap in the sorted values"),
            ], "5/6 agree with the paper; 1 of 2 TUNEL+ found"),
            ("Advanced", "55 min", [
                ("Import", "Bio-Formats; channel roles from emission metadata"),
                ("Segment", "StarDist trial rejected; scikit-image DoG "
                            "(s=4/30) + Otsu + watershed"),
                ("Filter", "area 200-5000 px, solidity >= 0.45, edges dropped"),
                ("Correct", "per-nucleus annulus 5-30 px, sigma-clipped mean"),
                ("Classify", "2-component GMM if justified, else median + 3*MAD"),
            ], "6/6 found, 0 invented; 5/6 agree"),
            ("Fast", "22 min", [
                ("Import", "Bio-Formats; overview ROI 608x608 at (565, 952)"),
                ("Segment", "StarDist 2D Versatile, prob 0.70, NMS 0.00"),
                ("Merge audit", "all 21 pairs tested, none merged -> 7 kept"),
                ("Correct", "per-nucleus annulus 6-20 px, DAPI-bright excluded"),
                ("Classify", "largest adjacent gap, conservative candidate"),
            ], "6/6 found + 1 invented; 5/6 agree"),
        ]),
    2: dict(
        title="Task 2 — proliferation index, single parasites in one frame of position 40",
        columns=[
            ("Published (manual)", "not recorded", [
                ("Draw ROIs", "by hand: each parasite + 3 local backgrounds, "
                              "4 subsets"),
                ("Measure", "ImageJ ROI Manager, red and green channels"),
                ("Export", "both ROI sets to Excel"),
                ("Compute", "corrected = raw - mean(3 bg); index = 10 - R/G"),
            ], "15 parasites, p = 0.0012"),
            ("Agent + human (repeat 1)", "not recorded", [
                ("Import", "Bio-Formats, series 4, final frame t = 262"),
                ("Pair the given ROIs", "9 parasites x 3 backgrounds per subset, "
                                        "1:1 verified"),
                ("Measure", "ImageJ, native 16-bit, red and green"),
                ("Deduplicate", "repeated controls collapsed: 36 rows -> 15 parasites"),
                ("Compute + test", "Python; index, Mann-Whitney"),
            ], "15/15 reproduced, p = 0.0012 — matches the paper"),
            ("Advanced", "85 min", [
                ("Import", "Bio-Formats: live series 4 + fixed series 39"),
                ("Register", "TurboReg / StackReg rigid, ~2-3 um residual"),
                ("Detect parasites", "micro_sam vit_b_lm, tiled 1024 px: "
                                     "105 proposals -> 5 kept"),
                ("Segment hosts", "Cellpose cyto3, auto diameter 72.6 px -> 6 hosts"),
                ("Classify hosts", "2-component GMM -> 2 TUNEL+, 4 TUNEL-"),
                ("Correct", "per-parasite annulus, inner 4-8, outer 10-11 px"),
            ], "3/15 found, 2 invented; no TUNEL+ parasite, so no test"),
            ("Fast", "69 min", [
                ("Import", "Bio-Formats"),
                ("Register", "TurboReg rigid on C1 brightfield / grid"),
                ("Detect parasites", "C2 + C3 evidence combined; 95 clumps "
                                     "excluded -> 22 kept"),
                ("Segment hosts", "Cellpose cyto3, auto diameter 35.9 px -> 7 hosts"),
                ("Assign", "no parasite fell inside any host mask"),
                ("Correct", "per-parasite annulus 4-12 px, widened to 18 px"),
            ], "0/15 found, 22 invented; every parasite Unassigned"),
        ]),
    3: dict(
        title="Task 3 — proliferation index, every parasite of every position",
        columns=[
            ("Published (manual)", "not recorded", [
                ("Draw ROIs", "by hand over 35 positions: 828 parasites, "
                              "3 backgrounds each"),
                ("Measure", "ImageJ ROI Manager, red and green, per position"),
                ("Export", "92 workbooks"),
                ("Compute", "index = 10 - R/G in Excel"),
            ], "828 parasites, p = 7.0e-05"),
            ("Agent + human (repeat 1)", "not recorded", [
                ("Inventory", "8 LIF files -> position, series and final-frame map"),
                ("Read the ROI sets", "top-level .roi only, nested RoiSet zips ignored"),
                ("Pair backgrounds", "one triplet per HOST, reused for its parasites"),
                ("Flag clumps", "ROI geometry + native-image evidence: 135 excluded"),
                ("Measure", "native 16-bit C2 / C3 per ROI"),
                ("Compute + test", "Mann-Whitney + host-median sensitivity"),
            ], "577/828 reproduced, p = 0.0078"),
            ("Advanced", "> 85 min", [
                ("Extract", "Bio-Formats: 8 LIFs -> 360 planes, 40 positions"),
                ("Register", "TurboReg rigid, 40/40 succeeded, rotation 1.72 deg"),
                ("Detect parasites", "Cellpose cyto3, diameter 24.3 px: "
                                     "0-5 per position"),
                ("STOPPED", "killed at the deadline before any measurement"),
            ], "no per-parasite table"),
            ("Fast", "> 85 min", [
                ("Extract", "Bio-Formats, 40 positions"),
                ("Detect parasites", "40/40 positions, ~43 400 candidates each "
                                     "-> 15 432 kept"),
                ("Register", "sampled 5 positions"),
                ("STOPPED", "killed at the deadline before any measurement"),
            ], "no per-parasite table"),
        ]),
}

STEP_WRAP, METHOD_WRAP = 26, 34
LINE = 0.0225                 # axes fraction per text line
PAD = 0.016


def draw_column(ax, x, width, name, time, steps, outcome):
    # "Published (manual)" is itself a key; "Agent + human (repeat 1)"
    # is the same configuration as "Agent + human" and shares its palette.
    style = MODES.get(name) or MODES[re.sub(r"\s*\(.*\)$", "", name)]
    ax.text(x, 0.965, name, ha="center", va="top", fontsize=12.5,
            fontweight="bold", color=style["edge"])
    # Time to result, as a pill under the name.
    ax.text(x, 0.918, time, ha="center", va="center", fontsize=10.5, color=INK,
            bbox=dict(boxstyle="round,pad=0.42", facecolor=style["fill"],
                      edgecolor=style["edge"], linewidth=1.1))

    y = 0.872
    for step, method in steps:
        s_lines = textwrap.wrap(step, STEP_WRAP) or [""]
        m_lines = textwrap.wrap(method, METHOD_WRAP)
        height = (len(s_lines) + len(m_lines)) * LINE + 2 * PAD
        stopped = step == "STOPPED"
        ax.add_patch(FancyBboxPatch(
            (x - width / 2, y - height), width, height,
            boxstyle="round,pad=0.004,rounding_size=0.012",
            facecolor="#f6f6f6" if stopped else style["fill"],
            edgecolor="#b0b0b0" if stopped else style["edge"],
            linewidth=1.3, linestyle=":" if stopped else "-",
            transform=ax.transData, clip_on=False, zorder=2))
        ty = y - PAD
        for line in s_lines:
            ax.text(x, ty, line, ha="center", va="top", fontsize=10.2,
                    fontweight="bold", color="#8a8a8a" if stopped else INK,
                    zorder=3)
            ty -= LINE
        for line in m_lines:
            ax.text(x, ty, line, ha="center", va="top", fontsize=8.8,
                    color="#6b6b6b", zorder=3)
            ty -= LINE
        y -= height
        if (step, method) != steps[-1]:
            ax.add_patch(FancyArrowPatch(
                (x, y), (x, y - 0.019), arrowstyle="-|>", mutation_scale=11,
                color="#9aa0a6", linewidth=1.1, clip_on=False, zorder=2))
        y -= 0.026



def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", default=".")
    p.add_argument("--tasks", nargs="+", type=int, default=[1, 2, 3])
    args = p.parse_args()

    for task in args.tasks:
        spec = TASKS[task]
        cols = spec["columns"]
        fig, ax = plt.subplots(figsize=(15.5, 9.4))
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.axis("off")
        width = 0.205
        centres = [(i + 0.5) / len(cols) for i in range(len(cols))]
        for x, (name, time, steps, outcome) in zip(centres, cols):
            draw_column(ax, x, width, name, time, steps, outcome)
        fig.suptitle(spec["title"], fontsize=14, y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.99])
        base = os.path.join(args.out_dir, f"figure_workflow_task{task}")
        fig.savefig(base + ".svg", facecolor="white")
        fig.savefig(base + ".png", dpi=200, facecolor="white")
        plt.close(fig)
        print(f"-> {base}.svg / .png   ({len(cols)} columns, "
              f"{max(len(c[2]) for c in cols)} steps deep)")


if __name__ == "__main__":
    main()
