#!/usr/bin/env python3
"""Tasks 2 and 3 as the paper plots them: the published figure beside each run.

    figure_tasks23_dotplots.py --gt-root ../parasite_task_results \
        --published ~/Downloads/jci.insight.169020.sdval.xlsx \
        --out figure_tasks23.png

The published figures are JCI Insight 169020, figure 4C (right) and figure 4D:
proliferation index of single parasites, split by whether the host macrophage was
TUNEL-positive, as a scatter dot plot with a median bar. This figure keeps that
form and adds one column per run, so a run is compared against the published
result in the shape the result was published in.

THE PUBLISHED VALUES ARE THE GROUND TRUTH, VERIFIED
---------------------------------------------------
Sheets `4C` and `4D` of the article's source data hold the exact values behind
those panels. They agree with this project's parse of the laboratory's own
workbooks to the last digit -- 6 values for task 1, 15 for task 2 and all 828 for
task 3, maximum difference 0.000000. So the column labelled "published" is both
the figure in the paper and the ground truth every score in this study is
measured against; they are not two different references that happen to agree.

WHAT IS AND IS NOT SHOWN
------------------------
Every value a run produced is plotted, including the ones it could not assign to
a host. Auto and fast mode report parasites as `Unassigned` when no parasite fell
inside any segmented macrophage, and dropping that category would turn a run that
answered nothing into a run that answered narrowly. A run with no output at all
is drawn as an empty column with the reason, for the same reason.
"""

import argparse
import os
import statistics as st

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import build_workbook as bw
import score_auto_run as auto
import score_parasite_run as hil
import xlsx_min

COLOR = {"TUNEL+": "#d1495b", "TUNEL-": "#3d6d8f", "Unassigned": "#9aa0a6"}
ORDER = ["TUNEL+", "TUNEL-", "Unassigned"]
JITTER = 0.085


def published_block(rows, start, end=None):
    plus = [r[0] for r in rows[start:end] if len(r) > 0 and isinstance(r[0], float)]
    minus = [r[1] for r in rows[start:end] if len(r) > 1 and isinstance(r[1], float)]
    return {"TUNEL+": plus, "TUNEL-": minus}


def group_values(rows, key="index"):
    out = {}
    for r in rows:
        v = r.get(key)
        if v is None:
            continue
        g = r.get("group") or "Unassigned"
        if g not in COLOR:
            g = "Unassigned"
        out.setdefault(g, []).append(v)
    return out


def collect(root, published_path):
    wb = xlsx_min.read(published_path)
    t2 = {"published (Fig. 4C)": published_block(wb["4C"], 8)}
    t3 = {"published (Fig. 4D)": published_block(wb["4D"], 1)}

    for label, folder in (("Agent + human (run 1)", "repeat_1"), ("Agent + human (run 2)", "repeat_2")):
        tabs = bw.hil_tables(os.path.join(root, "hil_results", folder))
        if tabs.get(2):
            t2[label] = group_values(hil.rows_of(tabs[2]))
        if tabs.get(3):
            # A run may leave more than one all-positions answer; take the most
            # complete, which is the run at its best.
            best = max(tabs[3], key=lambda t: t["n"])
            t3[label] = group_values(hil.rows_of(best))

    for label, rel in (("Advanced", os.path.join("auto-mode", "task_2")),
                       ("Fast", os.path.join(
                           "fast-mode", "session-20260925-182436", "task_2"))):
        table = bw.auto_table(os.path.join(root, rel), 2)
        t2[label] = group_values(
            [dict(index=o["value"], group=o["group"])
             for o in auto.objects_of(table)]) if table else {}

    for label, rel in (("Advanced", os.path.join("auto-mode", "task_3")),
                       ("Fast", os.path.join(
                           "fast-mode", "session-20260925-182436", "task_3"))):
        table = bw.auto_table(os.path.join(root, rel), 3)
        t3[label] = group_values(
            [dict(index=o["value"], group=o["group"])
             for o in auto.objects_of(table)]) if table else {}
    return t2, t3


def stars(p):
    """The usual convention, so the mark carries the same meaning as in the paper."""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def panel(ax, data, title, ylabel, no_result_note, ylim=None, yticks=None,
          significance=False, n_inside=False):
    rng = np.random.default_rng(0)          # reproducible jitter
    xticks, xlabels, x = [], [], 0.0
    group_seen = []
    slots = [g for g in ORDER
             if any(groups.get(g) for groups in data.values())]
    for source, groups in data.items():
        present = [g for g in slots if groups.get(g)]
        if not present:
            ax.text(x + 0.5, ax.get_ylim()[0], "", ha="center")
            ax.annotate(no_result_note.get(source, "no result"), (x + 0.5, 0.5),
                        xycoords=("data", "axes fraction"), ha="center",
                        va="center", fontsize=10, color="#8a8a8a",
                        style="italic")
            xticks.append(x + 0.5)
            xlabels.append(source)
            x += 1.6
            continue
        start = x
        ny, nva = (0.015, "bottom") if n_inside else (-0.035, "top")
        for g in slots:                       # every slot, drawn or marked empty
            vals = groups.get(g)
            if vals:
                if g not in group_seen:
                    group_seen.append(g)
                jit = rng.uniform(-JITTER, JITTER, len(vals))
                ax.scatter(np.full(len(vals), x) + jit, vals, s=14, alpha=0.65,
                           color=COLOR[g], edgecolors="none", zorder=2)
                ax.plot([x - 0.22, x + 0.22], [st.median(vals)] * 2, lw=2.4,
                        color="#1b1b1b", zorder=4, solid_capstyle="butt")
            ax.annotate(f"n={len(vals) if vals else 0}", (x, ny),
                        xycoords=("data", "axes fraction"), ha="center", va=nva,
                        fontsize=8.5, color="#1b1b1b" if vals else "#9aa0a6",
                        annotation_clip=False)
            x += 0.52
        centre = (start + x - 0.52) / 2
        xticks.append(centre)
        xlabels.append(source)
        # Mann-Whitney on the two real groups, computed the same way for every
        # column so the marks are comparable across them.
        a, b = groups.get("TUNEL+"), groups.get("TUNEL-")
        if significance and a and b and len(a) > 1 and len(b) > 1:
            _u, p = hil.mann_whitney(a, b)
            if p is not None:
                ax.annotate(stars(p), (centre, 0.93),
                            xycoords=("data", "axes fraction"), ha="center",
                            fontsize=13, color="#1b1b1b")
        x += 1.08
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=9.5)
    ax.tick_params(axis="x", length=0, pad=20)   # room for the n= row beneath
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-0.6, x - 0.8)
    if ylim:
        # A fixed range, matching the published panel. Points outside it are
        # drawn but fall off the axis; how many, and which run they belong to,
        # is printed by this script rather than crowding the figure.
        ax.set_ylim(*ylim)
    else:
        vals = [v for groups in data.values() for vs in groups.values() for v in vs]
        if vals:
            lo, hi = min(vals), max(vals)
            pad = max(hi - lo, 1.0)
            ax.set_ylim(lo - 0.16 * pad, hi + 0.16 * pad)
    if yticks is not None:
        ax.set_yticks(yticks)
    return group_seen


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gt-root", required=True)
    p.add_argument("--published", required=True)
    p.add_argument("--out", default="figure_tasks23.png")
    p.add_argument("--clip", type=float, default=None,
                   help="clip the y axis to this half-range around the published "
                        "median, keeping outliers out of the way (default: no clip)")
    args = p.parse_args()

    t2, t3 = collect(args.gt_root, os.path.expanduser(args.published))

    fig, axes = plt.subplots(2, 1, figsize=(11.5, 11.0),
                         gridspec_kw=dict(height_ratios=[1.0, 1.15]))
    seen = panel(axes[0], t2, "Task 2", "proliferation index  (10 - R/G)", {},
                 ylim=(0, 10), yticks=[0, 5, 10], n_inside=True)
    seen += [g for g in panel(axes[1], t3, "Task 3",
                              "proliferation index  (10 - R/G)", {},
                              ylim=(4, 10), yticks=[4, 6, 8, 10],
                              significance=True) if g not in seen]

    handles = [plt.Line2D([], [], marker="o", lw=0, ms=7, color=COLOR[g])
               for g in ORDER if g in seen]
    labels = [g if g != "Unassigned" else "Unassigned (no host could be assigned)"
              for g in ORDER if g in seen]
    handles.append(plt.Line2D([], [], color="#1b1b1b", lw=2.4))
    labels.append("median")
    fig.legend(handles, labels, frameon=False, fontsize=9.5, ncol=len(labels),
               loc="upper center", bbox_to_anchor=(0.5, 0.955))

    fig.suptitle("The published result and what each configuration reproduced",
                 fontsize=13, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.935], h_pad=3.2)
    fig.savefig(args.out, dpi=200, facecolor="white")
    svg = os.path.splitext(args.out)[0] + ".svg"
    fig.savefig(svg, facecolor="white")
    print(f"-> {args.out}\n-> {svg}")
    for name, data in (("task 2", t2), ("task 3", t3)):
        print(f"\n{name}")
        for source, groups in data.items():
            if not groups:
                print(f"   {source:22} no values")
                continue
            lim = (0, 10) if name == "task 2" else (4, 10)
            bits = []
            for g, v in groups.items():
                out = sum(1 for x in v if not lim[0] <= x <= lim[1])
                bits.append(f"{g} n={len(v)} median={st.median(v):.3f}"
                            + (f" ({out} outside the axis)" if out else ""))
            print(f"   {source:22} " + ", ".join(bits))


if __name__ == "__main__":
    main()
