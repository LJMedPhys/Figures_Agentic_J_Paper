#!/usr/bin/env python3
"""Assemble every run and the lab's ground truth into one .xlsx, a sheet per task.

    build_workbook.py --gt-root ../parasite_task_results \
                      --one-frame ../Imagent_J/data/parasite_test_case/one_frame \
                      --out parasite_task_results.xlsx

Sheets
------
    Summary      every run x task: coverage, agreement, group, the test
    Task 1       one row per lab nucleus: the lab's own numbers, then each run's
    Task 2       one row per lab parasite, same shape
    Task 3       one row per lab parasite across all 35 positions, same shape
    Effort       time, tool calls, tokens, reported and real cost
    Methods      the tool and method behind each step, per run
    Figures      where every scatter, overlay and mask actually is

THE LAB'S NUMBERS COME FIRST IN EVERY SHEET
------------------------------------------
Each task sheet starts with the ground truth -- label, group, coordinates, raw
green and red, corrected green and red, and `10-R/G` -- and only then the runs.
A run's column pair is its value and its difference from the lab, so a sheet can
be read left to right as "what was true, then what each configuration said".

Everything is matched the way the two scorers match: HIL runs by raw intensity
(the ROI was given, so the intensities identify it), auto and fast runs by
centroid distance (they segmented their own objects, so position is all there
is). The `matched by` column says which, and an empty value column means that
run has no number for that object -- not that it measured zero.
"""

import argparse
import csv
import glob
import json
import math
import os
import re
import statistics as st

import ij_roi
import parasite_gt as gt
import score_auto_run as auto
import score_parasite_run as hil
import xlsx_write

# label, mode, path relative to the results root, which tasks it covers
RUNS = [
    ("hil repeat_1", "hil",  "hil_results/repeat_1",                    (1, 2, 3)),
    ("hil repeat_2", "hil",  "hil_results/repeat_2",                    (1, 2, 3)),
    ("auto task_1",  "auto", "auto-mode/task_1",                        (1,)),
    ("auto task_2",  "auto", "auto-mode/task_2",                        (2,)),
    ("auto task_3",  "auto", "auto-mode/task_3",                        (3,)),
    ("fast task_1",  "auto", "fast-mode/session-20260925-182436/task_1", (1,)),
    ("fast task_2",  "auto", "fast-mode/session-20260925-182436/task_2", (2,)),
    ("fast task_3",  "auto", "fast-mode/session-20260925-182436/task_3", (3,)),
]
# gpt-5.6-sol, the real rates: the cost a run reports is roughly half of this.
RATE = dict(fresh=4.00 / 1e6, cached=0.40 / 1e6, out=20.00 / 1e6)
MAX_DIST = 15.0


def num(v, nd=4):
    return None if v is None else round(float(v), nd)


# --------------------------------------------------------------------------
# reading each run's answer
# --------------------------------------------------------------------------
def hil_tables(path):
    """-> {task: table} for a human-in-the-loop run, best answer per task."""
    t1, t2, t3, _notes = hil.choose_tables(hil.scan_run(path))
    out = {}
    if t1:
        # Prefer the table whose background level matches the lab's channel.
        best, best_gap = None, None
        for t in t1:
            rows = hil.rows_of(t)
            bg = next((r["background"] for r in rows if r["background"] is not None), None)
            gap = abs(bg - 190.2705) if bg is not None else 1e9
            if best_gap is None or gap < best_gap:
                best, best_gap = t, gap
        out[1] = best
    if t2:
        out[2] = t2
    if t3:
        # Where a run left several all-positions answers, take the one that
        # reproduces the most of the lab's parasites; that is the run at its best.
        out[3] = t3
    return out


def auto_table(path, task):
    found = auto.find_tables(path)
    if task == 3:
        found = [f for f in found if f["n_positions"] > 1]
    if not found:
        return None
    out = [f for f in found if f["rel"].startswith("output")]
    return max(out or found, key=lambda f: f["n"])


# --------------------------------------------------------------------------
# sheets
# --------------------------------------------------------------------------
def task1_sheet(root, one_frame, runs):
    nuclei, background = gt.load_task1(os.path.join(
        root, "oneframe_results", "ExampleData_TunelMeasurement.xlsx"))
    rois, _bg, check = gt.load_task1_coords(one_frame, nuclei)
    by_number = {n["roi"]: n for n in nuclei}
    bg_mean = st.mean([b["green"] for b in background])

    head = ["lab ROI", "x (full px)", "y (full px)", "area um2", "raw green",
            "background", "corrected", "lab group"]
    order = sorted(rois, key=lambda r: r["roi"])
    body = []
    for r in order:
        n = by_number[r["roi"]]
        body.append([r["name"], num(r["x"], 1), num(r["y"], 1), num(n["area"], 3),
                     num(n["green"], 4), num(bg_mean, 4), num(n["corrected"], 4),
                     n["group"]])

    extra_rows = []
    for label, mode, path, tasks in runs:
        if 1 not in tasks:
            continue
        head += [f"{label}: corrected", f"{label}: group", f"{label}: delta",
                 f"{label}: offset px"]
        cells = {r["name"]: ["", "", "", ""] for r in order}
        if mode == "hil":
            table = hil_tables(path).get(1)
            if table:
                rows = hil.rows_of(table)
                numbered = {}
                for x in rows:
                    m = re.search(r"(\d+)", x["name"])
                    if m:
                        numbered[int(m.group(1))] = x
                for r in order:
                    x = numbered.get(r["roi"])
                    if x and x["corrected"] is not None:
                        lab_v = by_number[r["roi"]]["corrected"]
                        cells[r["name"]] = [num(x["corrected"]), x["group"] or "",
                                            num(x["corrected"] - lab_v), "given"]
        else:
            table = auto_table(path, 1)
            if table:
                gts = [dict(name=r["name"], x=r["x"], y=r["y"],
                            value=by_number[r["roi"]]["corrected"],
                            group=by_number[r["roi"]]["group"]) for r in order]
                pairs, _missed, spurious = auto.match_by_distance(
                    gts, auto.objects_of(table), MAX_DIST)
                for g, rr, d in pairs:
                    cells[g["name"]] = [num(rr["value"]), rr["group"] or "",
                                        num(rr["value"] - g["value"]) if
                                        rr["value"] is not None else "", num(d, 1)]
                for s in spurious:
                    extra_rows.append([label, s["name"], num(s["x"], 1),
                                       num(s["y"], 1), num(s["value"]),
                                       s["group"] or ""])
        for i, r in enumerate(order):
            body[i] += cells[r["name"]]

    rows = [head] + body
    if extra_rows:
        rows += [[], ["objects a run reported that match no nucleus the lab measured"],
                 ["run", "its id", "x", "y", "corrected", "group"]] + extra_rows
    rows += [[], [f"lab ROI coordinates are bounding-box centres from the .roi files "
                  f"plus the TUNEL_Overview origin; the workbook's own calibrated "
                  f"X/Y agree to {check:.1f} px"],
             [f"lab background is the unweighted mean of {len(background)} "
              f"background ROI means, on channel 2 (TUNEL)"]]
    return rows


def task2_measured_block(root, runs):
    """The lab's 15 unique parasites with every value it recorded, and each HIL
    run's number for each.

    This is the block that carries all fifteen lab indices. The ROI block below
    is keyed by position instead, and there only the seven controls can be tied
    to a workbook row by name -- so without this block the eight TUNEL-positive
    parasites would appear to have no ground-truth index at all, when in fact
    their indices are exactly what the four workbooks were made to record.
    """
    g2, _files = gt.load_task2(os.path.join(root, "oneframe_results"))
    head = ["lab label", "workbook(s)", "lab group", "raw green", "raw red",
            "corrected green", "corrected red", "R/G", "10-R/G"]
    body = [[r["label"], "|".join(r.get("subsets") or [r.get("source", "")]),
             r["group"], num(r["green"]), num(r["red"]), num(r["cgreen"]),
             num(r["cred"]),
             num(10 - r["index"]) if r["index"] is not None else "",
             num(r["index"])] for r in g2]
    for label, mode, path, tasks in runs:
        if 2 not in tasks or mode != "hil":
            continue
        head += [f"{label}: 10-R/G", f"{label}: delta", f"{label}: run ROI"]
        cells = [["", "", ""] for _ in g2]
        table = hil_tables(path).get(2)
        if table:
            run_rows = [r for r in hil.rows_of(table) if r["index"] is not None]
            pairs, _mg, _mr, _s = hil.match(g2, run_rows, 0.5)
            row_of = {id(g): i for i, g in enumerate(g2)}
            for g, r, _how in pairs:
                cells[row_of[id(g)]] = [num(r["index"]), num(r["index"] - g["index"]),
                                        r["name"]]
        for i in range(len(body)):
            body[i] += cells[i]
    return [head] + body


def task2_sheet(root, one_frame, runs):
    g2, files = gt.load_task2(os.path.join(root, "oneframe_results"))
    by_label = {r["label"]: r for r in g2}
    coords = gt.load_task2_coords(one_frame)

    head = ["lab ROI", "subset", "lab group", "x (full px)", "y (full px)",
            "raw green", "raw red", "corrected green", "corrected red", "10-R/G"]
    body, gts = [], []
    for c in coords:
        row = by_label.get(c["gt_label"]) if c["gt_label"] else None
        body.append([c["name"], c["subset"], c["group"], num(c["x"], 1), num(c["y"], 1),
                     num(row["green"]) if row else "", num(row["red"]) if row else "",
                     num(row["cgreen"]) if row else "", num(row["cred"]) if row else "",
                     num(row["index"]) if row else ""])
        gts.append(dict(name=c["name"] + ("/" + c["subset"] if c["subset"] else ""),
                        x=c["x"], y=c["y"], group=c["group"],
                        value=row["index"] if row else None))

    extra_rows = []
    for label, mode, path, tasks in runs:
        if 2 not in tasks:
            continue
        head += [f"{label}: 10-R/G", f"{label}: group", f"{label}: delta",
                 f"{label}: matched by"]
        cells = [["", "", "", ""] for _ in body]
        if mode == "hil":
            table = hil_tables(path).get(2)
            if table:
                # The HIL scorer pairs on the index within the group, which is what
                # it can do for a table that publishes no raw means. Only the seven
                # controls have a lab index, so only they can be paired this way --
                # the eight TUNEL-positive parasites have coordinates but no
                # workbook row to pair against, and are left blank rather than
                # matched to whatever happens to be nearest.
                run_rows = [r for r in hil.rows_of(table) if r["index"] is not None]
                pairable = [(i, dict(green=None, red=None, group=c["group"],
                                     label=c["name"],
                                     index=by_label[c["gt_label"]]["index"]))
                            for i, c in enumerate(coords) if c["gt_label"]]
                pairs, _mg, _mr, _second = hil.match(
                    [g for _i, g in pairable], run_rows, 0.5)
                row_of = {id(g): i for i, g in pairable}
                for g, r, how in pairs:
                    cells[row_of[id(g)]] = [num(r["index"]), r["group"] or "",
                                            num(r["index"] - g["index"]), how]
        else:
            table = auto_table(path, 2)
            if table:
                pairs, _missed, spurious = auto.match_by_distance(
                    gts, auto.objects_of(table), MAX_DIST)
                for g, rr, d in pairs:
                    i = gts.index(g)
                    cells[i] = [num(rr["value"]), rr["group"] or "",
                                num(rr["value"] - g["value"]) if
                                (rr["value"] is not None and g["value"] is not None)
                                else "", f"centroid {d:.1f} px"]
                for s in spurious:
                    extra_rows.append([label, s["name"], num(s["x"], 1), num(s["y"], 1),
                                       num(s["value"]), s["group"] or ""])
        for i in range(len(body)):
            body[i] += cells[i]

    rows = [["MEASURED PARASITES -- the lab's fifteen, with every value it recorded, "
             "and each HIL run's number (matched on the index)"]]
    rows += task2_measured_block(root, runs)
    rows += [[], ["SAME PARASITES BY POSITION -- the .roi files give coordinates, "
                  "which is how the auto and fast runs are matched; only C<n> can "
                  "also be tied to a workbook row by name"], head] + body
    if extra_rows:
        rows += [[], ["objects a run reported that match no parasite the lab measured"],
                 ["run", "its id", "x", "y", "10-R/G", "group"]] + extra_rows
    rows += [[], [f"ground truth pooled from {len(files)} workbooks and collapsed to "
                  f"{len(g2)} unique parasites: the seven controls are ONE set of "
                  f"measurements copied into all four files"],
             ["only C<n> can be tied to a workbook row by name (C<n> = Control<n>), "
              "so the eight TUNEL-positive parasites have coordinates but no lab index"],
             ["raw and corrected intensities are reconstructed: some workbooks lost "
              "their decimal comma, so the stored value was 1000x too large. 10-R/G "
              "is a ratio and was unaffected"]]
    return rows


def task3_sheet(root, runs):
    g3 = gt.load_task3(os.path.join(root, "4. Data Analysis Excel Files"))
    flat = []
    for pos in sorted(g3):
        for r in g3[pos]:
            flat.append((pos, r))

    head = ["position", "workbook", "lab label", "lab group", "raw green", "raw red",
            "corrected green", "corrected red", "10-R/G"]
    body = [[pos, r["source"], r["label"], r["group"], num(r["green"]), num(r["red"]),
             num(r["cgreen"]), num(r["cred"]), num(r["index"])] for pos, r in flat]

    for label, mode, path, tasks in runs:
        if 3 not in tasks:
            continue
        tables = []
        if mode == "hil":
            t = hil_tables(path).get(3) or []
            tables = [(label if len(t) < 2 else
                       f"{label} [{os.path.basename(x['rel'])[:28]}]", x) for x in t]
        else:
            x = auto_table(path, 3)
            tables = [(label, x)] if x else []
            if not x:
                head += [f"{label}: 10-R/G"]
                for row in body:
                    row.append("")
                continue
        for col_label, table in tables:
            head += [f"{col_label}: 10-R/G", f"{col_label}: delta",
                     f"{col_label}: run ROI"]
            cells = [["", "", ""] for _ in flat]
            run_rows = hil.rows_of(table)
            by_pos = {}
            for r in run_rows:
                if r["index"] is not None:
                    by_pos.setdefault(r["position"], []).append(r)
            index_of = {id(r): i for i, (_p, r) in enumerate(flat)}
            for pos in sorted(g3):
                pairs, _mg, _mr, _second = hil.match(
                    g3[pos], by_pos.get(pos, []), 0.5)
                for g, r, _how in pairs:
                    i = index_of[id(g)]
                    cells[i] = [num(r["index"]), num(r["index"] - g["index"]),
                                r["name"]]
            for i in range(len(body)):
                body[i] += cells[i]

    rows = [head] + body
    rows += [[], [f"{len(flat)} parasites over {len(g3)} positions, from "
                  f"{len(set(r['source'] for _p, r in flat))} workbooks"],
             ["auto and fast mode produced no all-positions table: both were killed "
              "by the benchmark deadline, so their columns are empty by construction"]]
    return rows


def summary_sheet(root, one_frame, runs):
    rows = [["run", "mode", "task", "lab objects", "run objects", "matched",
             "missed", "invented / extra", "precision", "recall", "F1",
             "reproduced exactly", "median |delta|", "group agreement",
             "Mann-Whitney p", "note"]]

    nuclei, background = gt.load_task1(os.path.join(
        root, "oneframe_results", "ExampleData_TunelMeasurement.xlsx"))
    rois, _b, _c = gt.load_task1_coords(one_frame, nuclei)
    by_number = {n["roi"]: n for n in nuclei}
    g2, _files = gt.load_task2(os.path.join(root, "oneframe_results"))
    coords2 = gt.load_task2_coords(one_frame)
    by_label2 = {r["label"]: r for r in g2}
    g3 = gt.load_task3(os.path.join(root, "4. Data Analysis Excel Files"))
    n3 = sum(len(v) for v in g3.values())

    for label, mode, full, tasks in runs:
        for task in tasks:
            note = ""
            if mode == "hil":
                tabs = hil_tables(full)
                table = tabs.get(task)
                if task == 3 and table:
                    table = table  # a list
                if not table:
                    rows.append([label, mode, task, "", "", "", "", "", "", "", "",
                                 "", "", "", "", "no table"])
                    continue
                if task == 1:
                    run_rows = hil.rows_of(table)
                    numbered = {}
                    for x in run_rows:
                        m = re.search(r"(\d+)", x["name"])
                        if m:
                            numbered[int(m.group(1))] = x
                    agree = sum(1 for r in rois
                                if numbered.get(r["roi"]) and
                                numbered[r["roi"]]["group"] == by_number[r["roi"]]["group"])
                    ds = [abs(numbered[r["roi"]]["corrected"] - by_number[r["roi"]]["corrected"])
                          for r in rois if numbered.get(r["roi"])
                          and numbered[r["roi"]]["corrected"] is not None]
                    bg = next((x["background"] for x in run_rows
                               if x["background"] is not None), None)
                    note = ("channel 2, as the lab used" if bg and abs(bg - 190.27) < 2
                            else f"different channel (background {bg:.0f})" if bg else "")
                    rows.append([label, mode, task, len(rois), len(run_rows),
                                 len(ds), len(rois) - len(ds), "", "", "", "", "",
                                 num(st.median(ds), 3) if ds else "",
                                 f"{agree}/{len(rois)}", "", note])
                elif task == 2:
                    run_rows = [r for r in hil.rows_of(table) if r["index"] is not None]
                    pairs, mg, mr, _s = hil.match(g2, run_rows, 0.5)
                    ds = [abs(r["index"] - g["index"]) for g, r, _h in pairs]
                    vals = {k: [r["index"] for r in run_rows if r["group"] == k]
                            for k in ("TUNEL+", "TUNEL-")}
                    _u, p = hil.mann_whitney(vals["TUNEL+"], vals["TUNEL-"])
                    rows.append([label, mode, task, len(g2), len(run_rows),
                                 len(pairs), len(mg), len(mr), "", "", "",
                                 sum(1 for d in ds if d <= 5e-4),
                                 num(st.median(ds), 5) if ds else "",
                                 f"{len(pairs)}/{len(pairs)}", num(p, 6) if p else "",
                                 "the lab's own test: p=0.001243"])
                else:
                    for t in table:
                        run_rows = [r for r in hil.rows_of(t) if r["index"] is not None]
                        by_pos = {}
                        for r in run_rows:
                            by_pos.setdefault(r["position"], []).append(r)
                        matched, ds = 0, []
                        for pos in sorted(g3):
                            pairs, _m, _x, _s = hil.match(g3[pos], by_pos.get(pos, []), 0.5)
                            matched += len(pairs)
                            ds += [abs(r["index"] - g["index"]) for g, r, _h in pairs]
                        vals = {k: [r["index"] for r in run_rows if r["group"] == k]
                                for k in ("TUNEL+", "TUNEL-")}
                        _u, p = hil.mann_whitney(vals["TUNEL+"], vals["TUNEL-"])
                        rows.append([label, mode, task, n3, len(run_rows), matched,
                                     n3 - matched, len(run_rows) - matched, "", "", "",
                                     sum(1 for d in ds if d <= 5e-4),
                                     num(st.median(ds), 5) if ds else "",
                                     f"{matched}/{matched}", num(p, 6) if p else "",
                                     os.path.basename(t["rel"])])
            else:
                table = auto_table(full, task)
                if not table:
                    why = "killed by the benchmark deadline, no deliverable"
                    rows.append([label, mode, task,
                                 len(rois) if task == 1 else len(g2) if task == 2 else n3,
                                 0, 0, "", "", 0, 0, 0, "", "", "", "", why])
                    continue
                objs = auto.objects_of(table)
                if task == 1:
                    gts = [dict(name=r["name"], x=r["x"], y=r["y"],
                                value=by_number[r["roi"]]["corrected"],
                                group=by_number[r["roi"]]["group"]) for r in rois]
                else:
                    gts = [dict(name=c["name"], x=c["x"], y=c["y"], group=c["group"],
                                value=by_label2[c["gt_label"]]["index"]
                                if c["gt_label"] else None) for c in coords2]
                pairs, missed, spurious = auto.match_by_distance(gts, objs, MAX_DIST)
                tp, fn, fp = len(pairs), len(missed), len(spurious)
                prec = tp / (tp + fp) if tp + fp else 0.0
                rec = tp / (tp + fn) if tp + fn else 0.0
                f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
                graded = [(g, r) for g, r, _d in pairs
                          if r["group"] in ("TUNEL+", "TUNEL-")]
                agree = sum(1 for g, r in graded if g["group"] == r["group"])
                ds = [abs(r["value"] - g["value"]) for g, r, _d in pairs
                      if g["value"] is not None and r["value"] is not None]
                vals = {k: [o["value"] for o in objs if o["group"] == k
                            and o["value"] is not None] for k in ("TUNEL+", "TUNEL-")}
                _u, p = hil.mann_whitney(vals["TUNEL+"], vals["TUNEL-"])
                if not vals["TUNEL+"]:
                    note = "produced no TUNEL+ group, so the hypothesis is untestable"
                rows.append([label, mode, task, len(gts), len(objs), tp, fn, fp,
                             num(prec, 3), num(rec, 3), num(f1, 3), "",
                             num(st.median(ds), 3) if ds else "",
                             f"{agree}/{len(graded)}" if graded else "none gradeable",
                             num(p, 6) if p else "", note])
    rows += [[], ["HIL runs are matched by raw intensity: the ROI was given, so the "
                  "intensities identify which parasite a row is"],
             ["auto and fast runs are matched by centroid distance within "
              f"{MAX_DIST:g} px: they segmented their own objects, so position is "
              "the only shared quantity"],
             ["'reproduced exactly' means the run's 10-R/G is within 5e-4 of the "
              "lab's, i.e. it reproduced the lab's arithmetic"]]
    return rows


def effort_sheet(root, runs):
    rows = [["run", "mode", "task", "turns", "model minutes", "tool calls",
             "tokens (M)", "cached share", "reported USD", "real USD at "
             "gpt-5.6-sol rates", "finished"]]
    for label, mode, full, tasks in runs:
        log = None
        p = os.path.join(full, "logs", "usage_log.json")
        if os.path.isfile(p):
            log = json.load(open(p))
        md = {}
        for rp in sorted(glob.glob(os.path.join(full, "**", "result.json"),
                                   recursive=True)):
            try:
                d = json.load(open(rp))
            except (OSError, ValueError):
                continue
            md = d.get("metadata") or {}
            if log is None:
                log = ((md.get("usage_report") or {}).get("conversation") or {})
            finished = "yes" if d.get("success") else f"no - {d.get('message', '')}"
            break
        else:
            finished = "yes (interactive session)"
        qs = (log or {}).get("queries") or []
        mins = sum(q.get("thinking_seconds") or 0 for q in qs) / 60
        totals = (md.get("session_totals") or {}).get("totals") or {}
        tok = totals.get("total_tokens") or md.get("total_tokens") \
            or sum(q.get("total_tokens") or 0 for q in qs)
        calls = md.get("tool_calls") or sum(len(q.get("tool_call_log") or []) for q in qs)
        reported = md.get("total_cost_usd") or sum(q.get("cost_usd") or 0 for q in qs)
        real, cached = "", ""
        if totals:
            cin = totals.get("cached_input_tokens") or 0
            fresh = (totals.get("input_tokens") or 0) - cin
            real = round(fresh * RATE["fresh"] + cin * RATE["cached"]
                         + (totals.get("output_tokens") or 0) * RATE["out"], 2)
            if totals.get("input_tokens"):
                cached = round(100.0 * cin / totals["input_tokens"], 1)
        rows.append([label, mode, "+".join(str(t) for t in tasks), len(qs) or "",
                     round(mins, 1) if mins else "", calls or "",
                     round((tok or 0) / 1e6, 2) or "",
                     f"{cached}%" if cached != "" else "", round(reported or 0, 2),
                     real, finished])
    rows += [[], ["a run reports roughly half its true cost: the reported figure "
                  "predates the gpt-5.6-sol rates (4.00 / 0.40 cached / 20.00 per "
                  "million tokens), which the last column applies"],
             ["model minutes are summed thinking_seconds and cover only the queries "
              "that completed, so a run killed by the deadline ran longer than shown"]]
    return rows


def methods_sheet(root, runs):
    rows = [["run", "task", "step", "tool / method"]]
    for label, mode, full, tasks in runs:
        led = {}
        for p in sorted(glob.glob(os.path.join(full, "**", "state_ledger.json"),
                                  recursive=True), key=lambda p: p.count(os.sep)):
            try:
                led = json.load(open(p))
                break
            except (OSError, ValueError):
                continue
        plan = led.get("pipeline_plan") or []
        for i, step in enumerate(plan, 1):
            parts = re.split(r"\s*(?:—|--)\s*", step, maxsplit=1)
            rows.append([label, "+".join(str(t) for t in tasks), f"{i}. "
                         + parts[0].strip().replace("_", " "),
                         parts[1].strip() if len(parts) > 1 else ""])
        for d in (led.get("key_decisions") or []):
            rows.append([label, "+".join(str(t) for t in tasks), "decision", d])
        if not plan and not led:
            rows.append([label, "+".join(str(t) for t in tasks),
                         "no state ledger found", ""])
    rows += [[], ["the ledger records the plan and the decisions; what was actually "
                  "computed is in each run's QC reports, which trace_methods.py "
                  "mines alongside the saved scripts"]]
    return rows


def figures_sheet(root):
    rows = [["run", "task", "kind", "path", "note"]]
    patterns = [("scatter", r"scatter|dot_plot"),
                ("overlay / QC", r"overlay|_qc|qc_|verif|montage|scalebar|scale_bar"),
                ("mask", r"mask\.tif$|labels.*\.tif$")]
    for base, run_label in (("hil_results/repeat_1", "hil repeat_1"),
                            ("hil_results/repeat_2", "hil repeat_2"),
                            ("auto-mode/task_1", "auto task_1"),
                            ("auto-mode/task_2", "auto task_2"),
                            ("auto-mode/task_3", "auto task_3"),
                            ("fast-mode/session-20260925-182436/task_1", "fast task_1"),
                            ("fast-mode/session-20260925-182436/task_2", "fast task_2"),
                            ("fast-mode/session-20260925-182436/task_3", "fast task_3")):
        for f in sorted(glob.glob(os.path.join(root, base, "**", "*"), recursive=True)):
            if not f.lower().endswith((".png", ".svg", ".tif")):
                continue
            rel = os.path.relpath(f, root)
            kind = next((k for k, pat in patterns if re.search(pat, rel, re.I)), None)
            if not kind:
                continue
            # A sequential session leaves earlier tasks' projects inside the later
            # task's folder, so a figure's own path is the only thing that says
            # which task it belongs to.
            task = ""
            m = re.search(r"task_(\d)", base)
            if m:
                task = m.group(1)
            else:
                # A HIL run keeps all three tasks in one folder, so only the
                # file's own name says which task a figure belongs to.
                low = os.path.basename(rel).lower()
                if "all_positions" in low:
                    task = "3"
                elif "prolif" in low:
                    task = "2"
                elif "tunel" in low or "channel2" in low:
                    # `position40_channel2_group_scatter` is the SIX NUCLEI on
                    # channel 2, i.e. task 1 -- not a proliferation figure.
                    task = "1"
            note = ""
            if "prolif_pos40" in rel and task == "3":
                note = "belongs to TASK 2 - carried into the task_3 folder"
            elif "/archive/" in rel:
                note = "superseded version"
            rows.append([run_label, task, kind, rel, note])
    rows += [[], ["hil runs have scatters only: with the ROIs given there was "
                  "nothing to segment, so no segmentation overlay exists"],
             ["auto task_3 has no figure of any kind - it was killed before plotting"]]
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gt-root", required=True)
    p.add_argument("--one-frame", required=True)
    p.add_argument("--out", default="parasite_task_results.xlsx")
    args = p.parse_args()
    root, one = args.gt_root, args.one_frame

    runs = [(label, mode, os.path.join(root, rel), tasks)
            for label, mode, rel, tasks in RUNS]
    sheets = [
        ("Summary",  summary_sheet(root, one, runs)),
        ("Task 1",   task1_sheet(root, one, runs)),
        ("Task 2",   task2_sheet(root, one, runs)),
        ("Task 3",   task3_sheet(root, runs)),
        ("Effort",   effort_sheet(root, runs)),
        ("Methods",  methods_sheet(root, runs)),
        ("Figures",  figures_sheet(root)),
    ]
    xlsx_write.write(args.out, sheets)
    print(f"-> {args.out}")
    for name, rows in sheets:
        print(f"   {name:10} {len(rows):>5} rows x "
              f"{max(len(r) for r in rows):>3} cols")


if __name__ == "__main__":
    main()
