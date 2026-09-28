#!/usr/bin/env python3
"""Score an AUTO-MODE run: it had to find the objects, so this is a detection task.

    score_auto_run.py --run ../parasite_task_results/auto-mode/task_1 --task 1 \
                      --gt-root ../parasite_task_results \
                      --one-frame ../Imagent_J/data/parasite_test_case/one_frame

A DIFFERENT QUESTION FROM THE HIL SCORER
----------------------------------------
In the HIL condition the ROIs were given, so a run's raw intensities identify
which of the lab's parasites each row is, and `score_parasite_run.py` asks
whether the run reproduced the lab's number for it.

An auto-mode run segmented its own objects. Nothing in its output matches the
lab's numbers, and nothing is supposed to. So the questions become, in order:

  DETECTION  did it find the lab's objects at all, and how many did it invent?
             Matched by centroid distance, because position is the only thing the
             two analyses share. This is reported first because everything else
             is conditional on it -- an index is meaningless for an object that
             is not the one the lab measured.
  GROUP      of the objects it did find, did it put them in the same group?
  VALUE      and did it get a comparable number? Only where the lab's own label
             for that object is unambiguous: `C<n>.roi` is `Control<n>` in the
             workbooks, so those seven have a ground-truth index. The eight
             TUNEL-positive parasites cannot be tied to a workbook row by name,
             so their values are compared as a distribution, not one by one.

The lab's coordinates come from the .roi files, read directly (see ij_roi); the
workbooks hold no coordinates. Task 1's ROIs are stored relative to
`TUNEL_Overview.roi` and the origin is added; the check that this is right is
that the workbook's own calibrated X/Y then agree to ~2 px.
"""

import argparse
import csv
import math
import os
import re
import statistics as st

import ij_roi
import parasite_gt as gt

CENTROID_X = ("centroid_x", "x", "center_x", "centre_x", "cx", "x_px")
CENTROID_Y = ("centroid_y", "y", "center_y", "centre_y", "cy", "y_px")
VALUE_COLUMNS = ("proliferation_index", "corrected_tunel_intensity",
                 "corrected_tunel_mean", "corrected_tunel_signal",
                 "corrected_intensity", "index")
GROUP_COLUMNS = ("group", "group_label", "tunel_group", "class", "tunel_class",
                 "host_class")
SKIP_DIRS = {"scripts", "raw_images", "processed_images", "figures", "logs"}


def _pick(fields, candidates):
    low = {(f or "").strip().lower(): f for f in fields or []}
    for c in candidates:
        if c in low:
            return low[c]
    return None


def _float(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def _norm_group(v):
    s = str(v or "").strip().upper().replace(" ", "").replace("_", "")
    if s in ("TUNEL+", "TUNELPOS", "TUNELPOSITIVE", "POSITIVE", "POS", "APOPTOTIC"):
        return "TUNEL+"
    if s in ("TUNEL-", "TUNELNEG", "TUNELNEGATIVE", "NEGATIVE", "NEG", "CONTROL",
             "NONAPOPTOTIC"):
        return "TUNEL-"
    return str(v or "").strip() or None


def find_tables(run_dir):
    """Every CSV with per-object centroids and a value — that is the answer table."""
    found = []
    for dirpath, dirnames, filenames in os.walk(run_dir):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in sorted(filenames):
            if not fn.lower().endswith(".csv"):
                continue
            path = os.path.join(dirpath, fn)
            try:
                with open(path, newline="") as fh:
                    rd = csv.DictReader(fh)
                    rows = list(rd)
                    fields = rd.fieldnames
            except (OSError, UnicodeDecodeError, csv.Error):
                continue
            xc, yc = _pick(fields, CENTROID_X), _pick(fields, CENTROID_Y)
            vc = _pick(fields, VALUE_COLUMNS)
            if not (rows and xc and yc and vc):
                continue
            pos_col = _pick(fields, ("position", "global_position", "position_tag"))
            npos = len({re.search(r"(\d+)", str(r.get(pos_col) or "")).group(1)
                        for r in rows
                        if pos_col and re.search(r"(\d+)", str(r.get(pos_col) or ""))}) \
                if pos_col else 0
            found.append(dict(path=path, rel=os.path.relpath(path, run_dir),
                              rows=rows, n=len(rows), x=xc, y=yc, value=vc,
                              position_col=pos_col, n_positions=npos,
                              group=_pick(fields, GROUP_COLUMNS),
                              name=_pick(fields, ("parasite_id", "nucleus_id",
                                                  "id", "label", "roi_name"))))
    return found


def objects_of(table):
    out = []
    for r in table["rows"]:
        x, y = _float(r.get(table["x"])), _float(r.get(table["y"]))
        if x is None or y is None:
            continue
        out.append(dict(x=x, y=y, value=_float(r.get(table["value"])),
                        group=_norm_group(r.get(table["group"])) if table["group"] else None,
                        name=str(r.get(table["name"]) or "").strip()))
    return out


def match_by_distance(gts, runs, max_dist):
    """Optimal 1:1 assignment on centroid distance, rejecting pairs beyond the cap."""
    if not gts or not runs:
        return [], list(gts), list(runs)
    cost = [[math.hypot(g["x"] - r["x"], g["y"] - r["y"]) for r in runs] for g in gts]
    try:
        import numpy as np
        from scipy.optimize import linear_sum_assignment
        ri, ci = linear_sum_assignment(np.asarray(cost))
        pairs = [(int(i), int(j), cost[i][j]) for i, j in zip(ri, ci)]
    except ImportError:
        cand = sorted((c, i, j) for i, row in enumerate(cost) for j, c in enumerate(row))
        gi, gj, pairs = set(), set(), []
        for c, i, j in cand:
            if i in gi or j in gj:
                continue
            gi.add(i), gj.add(j)
            pairs.append((i, j, c))
    kept = [(gts[i], runs[j], d) for i, j, d in pairs if d <= max_dist]
    mg = {id(g) for g, _r, _d in kept}
    mr = {id(r) for _g, r, _d in kept}
    return (kept, [g for g in gts if id(g) not in mg],
            [r for r in runs if id(r) not in mr])


def describe(values):
    if not values:
        return "n=0"
    s = f"n={len(values)} median={st.median(values):.3f} mean={st.mean(values):.3f}"
    if len(values) > 1:
        s += f" sd={st.stdev(values):.3f}"
    return s + f" range=[{min(values):.3f}, {max(values):.3f}]"


def report(title, gts, runs, max_dist, value_label):
    print("=" * 78)
    print(title)
    print("=" * 78)
    pairs, missed, spurious = match_by_distance(gts, runs, max_dist)
    tp, fn, fp = len(pairs), len(missed), len(spurious)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    print(f"  the lab measured {len(gts)} object(s); the run reported {len(runs)}")
    print(f"\n  DETECTION at {max_dist:g} px")
    print(f"     found      {tp}")
    print(f"     missed     {fn}")
    print(f"     invented   {fp}")
    print(f"     precision {prec:.3f}   recall {rec:.3f}   F1 {f1:.3f}")
    if pairs:
        d = sorted(p[2] for p in pairs)
        print(f"     centroid offset of matches: median {st.median(d):.1f} px, "
              f"max {max(d):.1f} px")
    for g in missed:
        print(f"     MISSED   {g['name']:22} at ({g['x']:.0f}, {g['y']:.0f})"
              + (f"  lab {value_label} {g['value']:.4f}" if g.get("value") is not None else ""))
    for r in spurious:
        print(f"     INVENTED {str(r['name'])[:22]:22} at ({r['x']:.0f}, {r['y']:.0f})"
              + (f"  run {value_label} {r['value']:.4f}" if r.get("value") is not None else "")
              + (f"  [{r['group']}]" if r.get("group") else ""))

    graded = [(g, r, d) for g, r, d in pairs
              if g.get("group") and r.get("group") in ("TUNEL+", "TUNEL-")]
    if graded:
        agree = sum(1 for g, r, _d in graded if g["group"] == r["group"])
        print(f"\n  GROUP      {agree}/{len(graded)} of the objects it found are in the "
              f"same group as the lab put them")
        for g, r, _d in graded:
            if g["group"] != r["group"]:
                print(f"     {g['name']:22} lab {g['group']:7} run {r['group']}")
    ungraded = [r for _g, r, _d in pairs if r.get("group") not in ("TUNEL+", "TUNEL-")]
    if ungraded:
        labels = sorted({str(r.get("group")) for r in ungraded})
        print(f"     ({len(ungraded)} matched object(s) carry no usable group: "
              f"{', '.join(labels)})")

    valued = [(g, r) for g, r, _d in pairs
              if g.get("value") is not None and r.get("value") is not None]
    if valued:
        deltas = [r["value"] - g["value"] for g, r in valued]
        print(f"\n  VALUE      {value_label} on the {len(valued)} object(s) whose lab "
              f"row is identifiable")
        print(f"     median |delta| {st.median([abs(x) for x in deltas]):.3f}, "
              f"max {max(abs(x) for x in deltas):.3f}")
        for (g, r), d in zip(valued, deltas):
            print(f"     {g['name']:22} lab {g['value']:8.4f}  run {r['value']:8.4f}  "
                  f"{d:+8.4f}")

    print(f"\n  CONCLUSION")
    for src, data in (("lab", gts), ("run", runs)):
        for grp in ("TUNEL+", "TUNEL-"):
            vals = [o["value"] for o in data
                    if o.get("group") == grp and o.get("value") is not None]
            if vals:
                print(f"     {src:4} {grp:7} {describe(vals)}")
    print()
    return dict(tp=tp, fn=fn, fp=fp, precision=prec, recall=rec, f1=f1)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True)
    p.add_argument("--task", required=True, type=int, choices=[1, 2, 3])
    p.add_argument("--gt-root", required=True)
    p.add_argument("--one-frame", required=True,
                   help="the dataset's one_frame/ folder, which holds the ROI files")
    p.add_argument("--max-dist", type=float, default=15.0,
                   help="centroid distance beyond which a pair is not the same "
                        "object (default 15 px; a parasite is about 10 px across, "
                        "a nucleus about 30)")
    p.add_argument("--table", help="score this CSV instead of the discovered one")
    args = p.parse_args()

    print("=" * 78)
    print(f"AUTO-MODE RUN  {args.run}   task {args.task}")
    print("=" * 78)
    found = find_tables(args.run)
    print(f"  {len(found)} table(s) with per-object centroids and a value:")
    for f in found:
        print(f"      {f['rel']:58.58} {f['n']:>4} rows   "
              f"value={f['value']}")
    # A run that did several tasks in one session leaves the earlier tasks'
    # projects in the later task's folder. Scoring task 3 against a single-position
    # table carried over from task 2 would silently report task 2's answer as
    # task 3's, so task 3 requires a table that actually spans positions.
    if args.task == 3:
        multi = [f for f in found if f["n_positions"] > 1]
        if not multi:
            carried = [f for f in found if f["n_positions"] <= 1]
            if carried:
                print(f"\n  !! the {len(carried)} table(s) above cover a single "
                      f"position. They belong to an earlier task whose project "
                      f"folder was carried into this one, not to task 3.")
        found = multi
    if args.table:
        table = next((f for f in found if f["path"].endswith(args.table)
                      or f["rel"] == args.table), None)
        if table is None:
            raise SystemExit(f"--table {args.table} is not one of the tables found")
    else:
        # The deliverable is the largest per-object table under output/ if there is
        # one, else the largest anywhere: a run's `output/` folder is what the task
        # asked it to produce, and its data/ folder holds the working tables.
        out = [f for f in found if f["rel"].startswith("output")]
        table = max(out or found, key=lambda f: f["n"], default=None)
    if table is None:
        print("\n  !! no per-object table with centroids exists in this run.")
        print("     Nothing can be scored: the run did not produce the deliverable.")
        return
    print(f"\n  scoring: {table['rel']}")
    runs = objects_of(table)
    print()

    if args.task == 1:
        nuclei, _bg = gt.load_task1(os.path.join(
            args.gt_root, "oneframe_results", "ExampleData_TunelMeasurement.xlsx"))
        rois, _bgr, check = gt.load_task1_coords(args.one_frame, nuclei)
        by_number = {n["roi"]: n for n in nuclei}
        gts = [dict(name=r["name"], x=r["x"], y=r["y"],
                    value=by_number[r["roi"]]["corrected"],
                    group=by_number[r["roi"]]["group"]) for r in rois]
        print(f"  ground truth: {len(gts)} nuclei, located from the ROI files; "
              f"the workbook's own X/Y agree to {check:.1f} px")
        report("TASK 1 -- did it find and classify the lab's six nuclei?",
               gts, runs, args.max_dist, "corrected TUNEL")

    elif args.task == 2:
        g2, _files = gt.load_task2(os.path.join(args.gt_root, "oneframe_results"))
        by_label = {r["label"]: r for r in g2}
        coords = gt.load_task2_coords(args.one_frame)
        gts = []
        for c in coords:
            row = by_label.get(c["gt_label"]) if c["gt_label"] else None
            gts.append(dict(name=c["name"] + (f"/{c['subset']}" if c["subset"] else ""),
                            x=c["x"], y=c["y"], group=c["group"],
                            value=row["index"] if row else None))
        known = sum(1 for g in gts if g["value"] is not None)
        print(f"  ground truth: {len(gts)} unique parasites located from the ROI "
              f"files; {known} of them ({', '.join(sorted(by_label))}) can be tied "
              f"to a workbook row by name, so only those have a lab index")
        report("TASK 2 -- did it find the lab's parasites in the position-40 frame?",
               gts, runs, args.max_dist, "proliferation index")

    else:
        print("  task 3 spans all 40 positions; its ground truth is the per-position "
              "ROI sets and workbooks. Scoring it needs a per-parasite table with "
              "positions, which this run must contain for the comparison to mean "
              "anything.")
        report("TASK 3 -- all positions", [], runs, args.max_dist,
               "proliferation index")


if __name__ == "__main__":
    main()
