#!/usr/bin/env python3
"""The task-1 comparison figure: what each condition drew, on the same image.

    figure_task1_overlays.py --gt-root ../parasite_task_results \
                             --one-frame ../Imagent_J/data/parasite_test_case/one_frame \
                             --out figure_task1.png

WHY THIS IS THE FIGURE WORTH SHOWING
------------------------------------
Task 1 is the only task where all four conditions produced a comparable answer,
and the only one where the human-in-the-loop advantage is nil: the six nuclei are
given ovals either way. So it is the one place where the conditions can be put
side by side on the same pixels and the difference read off directly.

Top row, all on the same DAPI crop:
  the lab's six ROIs, then the objects auto mode segmented, then fast mode's.
  Outline colour is the TUNEL call, so a disagreement is visible as a colour
  change on the same nucleus rather than as a number in a table.

Bottom row:
  every nucleus's corrected TUNEL signal, lab against run, each point ringed in
  red where its own source called it TUNEL+. This is where the figure earns its
  place: the three traces track each other, so the runs agree with the lab on the
  MEASUREMENT, while Cell 6 carries a ring for the lab alone -- they disagree only
  on where to cut. No threshold derived from these six values recovers the lab's
  second positive, its own included.

The masks are full-frame label images; the lab's ROIs are ovals read from the
.roi files, offset by the overview origin. Nothing here is recomputed -- every
number comes from the run's own output.
"""

import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf
from matplotlib.patches import Ellipse
from skimage.measure import find_contours

import ij_roi
import parasite_gt as gt

POS = dict(color_plus="#d1495b", color_minus="#3d6d8f", edge="#1b1b1b")


def stretch(img, lo=1.0, hi=99.7):
    """Percentile contrast stretch — these are 16-bit planes with a narrow range."""
    a, b = np.percentile(img, [lo, hi])
    return np.clip((img.astype(np.float32) - a) / max(b - a, 1e-6), 0, 1)


def load_run(path):
    """-> {id: (x, y, corrected, group)} from a run's per-nucleus CSV."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            key = next((k for k in r if k.lower() in
                        ("nucleus_id", "id", "roi_name", "label")), None)
            out[str(r[key]).strip()] = (
                float(r["centroid_x"]), float(r["centroid_y"]),
                float(r["corrected_tunel_intensity"]), r["group"].strip())
    return out


def outline(ax, mask, ids_groups, origin, colors):
    """Draw each label's contour in the colour of its TUNEL call."""
    ox, oy = origin
    for lid, group in ids_groups.items():
        m = mask == lid
        if not m.any():
            continue
        for c in find_contours(m.astype(float), 0.5):
            ax.plot(c[:, 1] - ox, c[:, 0] - oy, lw=1.6,
                    color=colors["color_plus"] if group == "TUNEL+"
                    else colors["color_minus"])


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gt-root", required=True)
    p.add_argument("--one-frame", required=True)
    p.add_argument("--out", default="figure_task1.png")
    args = p.parse_args()
    root = args.gt_root

    nuclei, background = gt.load_task1(os.path.join(
        root, "oneframe_results", "ExampleData_TunelMeasurement.xlsx"))
    rois, _bg, _check = gt.load_task1_coords(args.one_frame, nuclei)
    by_number = {n["roi"]: n for n in nuclei}
    ov = ij_roi.read(os.path.join(args.one_frame, "TUNEL_Overview.roi"))
    ox, oy, w, h = ov["left"], ov["top"], ov["w"], ov["h"]

    base = os.path.join(root, "auto-mode", "task_1",
                        "leishmania_tunel_position40_task1", "raw_images")
    dapi = tf.imread(os.path.join(base, "position40_DAPI_raw16bit_ch3.tif"))
    tunel = tf.imread(os.path.join(base, "position40_TUNEL_raw16bit_ch2.tif"))
    dapi_c = stretch(dapi[oy:oy + h, ox:ox + w])
    tunel_c = stretch(tunel[oy:oy + h, ox:ox + w], 2.0, 99.9)

    auto_csv = os.path.join(root, "auto-mode", "task_1", "tunel_nuclei.csv")
    fast_csv = os.path.join(root, "fast-mode", "session-20260925-182436",
                            "task_1", "tunel_nuclei.csv")
    auto = load_run(auto_csv)
    fast = load_run(fast_csv)
    auto_mask = tf.imread(os.path.join(root, "auto-mode", "task_1",
                                       "tunel_nuclei_mask.tif"))
    fast_mask = tf.imread(os.path.join(root, "fast-mode",
                                       "session-20260925-182436", "task_1",
                                       "tunel_nuclei_mask.tif"))

    fig = plt.figure(figsize=(13.5, 8.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.22, wspace=0.06)

    panels = [
        ("Published: 6 ROIs drawn by hand", dapi_c, "roi"),
        (f"Advanced: {len(auto)} objects segmented", dapi_c, "auto"),
        (f"Fast: {len(fast)} objects segmented", dapi_c, "fast"),
    ]
    for k, (title, img, kind) in enumerate(panels):
        ax = fig.add_subplot(gs[0, k])
        ax.imshow(img, cmap="gray", interpolation="nearest")
        if kind == "roi":
            for r in sorted(rois, key=lambda r: r["roi"]):
                n = by_number[r["roi"]]
                ax.add_patch(Ellipse((r["x"] - ox, r["y"] - oy), r["w"], r["h"],
                                     fill=False, lw=1.8,
                                     edgecolor=POS["color_plus"] if n["group"] == "TUNEL+"
                                     else POS["color_minus"]))
                ax.annotate(str(r["roi"]), (r["x"] - ox + 18, r["y"] - oy - 14),
                            color="w", fontsize=8)
        else:
            mask, runs = (auto_mask, auto) if kind == "auto" else (fast_mask, fast)
            outline(ax, mask, {int(k2): v[3] for k2, v in runs.items()
                               if k2.isdigit()}, (ox, oy), POS)
            for k2, (x, y, _v, _g) in runs.items():
                ax.annotate(k2, (x - ox + 18, y - oy - 14), color="w", fontsize=8)
        ax.set_title(title, fontsize=10.5)
        ax.set_xlim(0, w); ax.set_ylim(h, 0)
        ax.set_xticks([]); ax.set_yticks([])
    # A scale bar on the first panel only; 0.32231 um per pixel.
    ax0 = fig.axes[0]
    bar = 20 / 0.32231030776746455
    ax0.plot([20, 20 + bar], [h - 24, h - 24], lw=3, color="w")
    ax0.annotate("20 um", (20, h - 34), color="w", fontsize=8)

    # ---- bottom: the measurement, and where each run cut it -----------------
    ax = fig.add_subplot(gs[1, :])
    order = sorted(rois, key=lambda r: r["roi"])
    xs = np.arange(len(order))
    lab_v = [by_number[r["roi"]]["corrected"] for r in order]
    lab_g = [by_number[r["roi"]]["group"] for r in order]

    def series(runs, mask):
        """Each run's value for the lab's nucleus n, matched by nearest centroid."""
        vals, grps = [], []
        for r in order:
            best, bd = None, 1e9
            for _k, (x, y, v, g) in runs.items():
                d = (x - r["x"]) ** 2 + (y - r["y"]) ** 2
                if d < bd:
                    best, bd = (v, g), d
            vals.append(best[0] if bd <= 15 ** 2 else np.nan)
            grps.append(best[1] if bd <= 15 ** 2 else "")
        return vals, grps

    auto_v, auto_g = series(auto, auto_mask)
    fast_v, fast_g = series(fast, fast_mask)

    ax.plot(xs, lab_v, "o-", color="#1b1b1b", lw=1.4, ms=9, label="Published", zorder=3)
    ax.plot(xs, auto_v, "s--", color="#2a9d8f", lw=1.2, ms=8, label="Advanced")
    ax.plot(xs, fast_v, "^--", color="#e07a1f", lw=1.2, ms=8, label="Fast")

    # A red ring on every point its own source called TUNEL+, using the same
    # colour as the TUNEL+ outlines in the panels above. This carries the
    # classification on the marker itself, so the one disagreement -- Cell 6,
    # where only the lab has a ring -- is read off the figure rather than from a
    # caption, and the same mark means the same thing in both rows.
    for values, groups, marker in ((lab_v, lab_g, "o"), (auto_v, auto_g, "s"),
                                   (fast_v, fast_g, "^")):
        hit = [(x, v) for x, v, g in zip(xs, values, groups) if g == "TUNEL+"]
        if hit:
            ax.scatter([x for x, _v in hit], [v for _x, v in hit], marker=marker,
                       s=135, facecolors="none", edgecolors=POS["color_plus"],
                       linewidths=1.8, zorder=5)
    ax.set_xlim(-0.45, len(order) - 0.55)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"Cell {r['roi']}" for r in order])
    ax.set_ylabel("corrected TUNEL signal (a.u.)")
    ax.set_title("Tunel singal intensities for different modes", fontsize=10.5)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(plt.Line2D([], [], marker="o", ms=11, lw=0, markerfacecolor="none",
                              markeredgecolor=POS["color_plus"], markeredgewidth=1.8))
    labels.append("called TUNEL+ by that source")
    ax.legend(handles, labels, frameon=False, fontsize=9, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)

    #fig.suptitle("Task 1 — TUNEL classification of six nuclei, position 40",
                 #fontsize=13, y=0.985)
    fig.savefig(args.out, dpi=200, bbox_inches="tight", facecolor="white")
    svg = os.path.splitext(args.out)[0] + ".svg"
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    print(f"-> {args.out}\n-> {svg}")
    print(f"   lab  {['%.2f' % v for v in lab_v]}  {lab_g}")
    print(f"   auto {['%.2f' % v for v in auto_v]}  {auto_g}")
    print(f"   fast {['%.2f' % v for v in fast_v]}  {fast_g}")


if __name__ == "__main__":
    main()
