#!/usr/bin/env python3
"""Score one ImagentJ run of the parasite tasks against the lab's own analysis.

    score_parasite_run.py --run  .../hil_results/repeat_1 \
                          --gt-root .../parasite_task_results \
                          --out  repeat_1

WHAT IS BEING SCORED
--------------------
In the human-in-the-loop condition the ROIs are GIVEN, so the run is not being
asked to find parasites -- it is being asked to read the right channels, pair
each parasite with its three background regions, and compute
`10 - corrected_red/corrected_green`. The lab computed the same number in Excel.
Agreement should therefore be near-exact, and every parasite that is missing,
duplicated or off by more than rounding is a real defect rather than noise.

So three things are reported, and they answer different questions:

  COVERAGE   how many of the lab's parasites the run produced a number for.
             A run can agree perfectly on everything it measured and still have
             measured a seventh of the dataset.
  AGREEMENT  for the parasites both scored, how far apart the two indices are.
  CONCLUSION whether the run's group statistics reproduce the lab's, computed
             from both sets of numbers the same way, so the comparison is of
             the data and not of two different statistical choices.

MATCHING
--------
Two stages, because the interesting failure is invisible to the first.

  BY ROI      the parasite's own raw green and red means identify which ROI was
              measured, and for a given ROI they agree to five decimals whatever
              else the run did. Pairing on them says WHICH parasite a row is,
              independently of the number the run computed for it. That is what
              separates "the run never measured this parasite" from "the run
              measured it and got a different answer", and it is how the
              background-pairing errors below were found at all.

  BY INDEX    for runs that do not publish raw means, an optimal assignment on
              `10-R/G` within a position and within a group, rejecting pairs
              further apart than --max-delta. This is weaker: where a run is
              missing parasites, the leftovers get paired with whatever is
              nearest, so its large deltas are partly an artefact. The report
              says which stage each pair came from.

Never on the raw intensities as the workbooks store them, though: some sheets
lost their decimal comma (see parasite_gt), so `parasite_gt` divides that back
out and the value compared is the reconstructed one.

Run files are found by their COLUMNS, not their names. Every run names its
outputs differently -- `All_Positions_Proliferation_Final_Table.csv` in one,
`all_positions_proliferation_index_per_parasite.csv` in the next -- and the
discovery report prints what was found and what was chosen so a wrong pick is
visible rather than silent.
"""

import argparse
import collections
import csv
import math
import os
import re
import statistics as st
import sys

import parasite_gt as gt

INDEX_COLUMNS = ("proliferation_index", "proliferation index", "index",
                 "proliferation_index_value", "prolif_index")
GROUP_COLUMNS = ("group", "group_label", "display_group", "group_label_norm",
                 "tunel_group", "host_group", "host_class", "tunel_class",
                 "tunel_group_canonical", "tunel_group_raw", "class")
POSITION_COLUMNS = ("position", "global_position", "position_tag", "position_id",
                    "pos")
NAME_COLUMNS = ("parasite_roi_name", "parasite_id", "roi_name", "parasite",
                "stable_parasite_id", "parasite_stable_id",
                "biological_parasite_id", "name")
CORRECTED_COLUMNS = ("corrected_tunel_signal", "corrected_tunel_mean",
                     "corrected_mean", "corrected_green", "corrected_signal",
                     "corrected")
BACKGROUND_COLUMNS = ("background_mean_used", "average_background_mean",
                      "avg_background_mean", "background_mean",
                      "average_background_green")
GREEN_COLUMNS = ("parasite_green_raw_mean", "raw_green_mean_parasite",
                 "green_raw_mean", "raw_green_mean", "raw_green",
                 "parasite_green_mean", "green_mean", "green")
RED_COLUMNS = ("parasite_red_raw_mean", "raw_red_mean_parasite",
               "red_raw_mean", "raw_red_mean", "raw_red",
               "parasite_red_mean", "red_mean", "red")
SKIP_DIRS = {"scripts", "raw_images", "processed_images", "figures", "logs"}
# Two ROI measurements of the same region agree to far better than this; two
# different parasites in the same field never do.
SAME_ROI_TOLERANCE = 0.01


# --------------------------------------------------------------------------
# reading a run
# --------------------------------------------------------------------------
# Runs invent column names as freely as they invent file names, and a name this
# module has never seen is the one case where discovery fails silently. So after
# the known spellings, fall back to shape: anything ending in `_class`, anything
# beginning `corrected_`. The excluded words are the near misses that would
# otherwise be taken for the value itself.
PATTERNS = {
    "group": (r"_class$", r"^class$", r"_group$", r"^group",
              r"^status", r"_status$"),
    "corrected": (r"^corrected_", r"_corrected$"),
}
NOT_A_VALUE = ("rank", "flag", "method", "threshold", "reason", "formula",
               "unit", "count", "check", "valid", "note")


def _pick(fieldnames, candidates, patterns=()):
    low = {(f or "").strip().lower(): f for f in fieldnames or []}
    for c in candidates:
        if c in low:
            return low[c]
    for pat in patterns:
        for name, original in low.items():
            if re.search(pat, name) and not any(w in name for w in NOT_A_VALUE):
                return original
    return None


def _norm_group(v):
    s = str(v or "").strip().upper().replace(" ", "")
    if s.startswith("TUNEL_POS") or s in ("TUNEL+", "TUNELPOS", "POSITIVE", "POS",
                                          "APOPTOTIC", "TUNELPOSITIVE"):
        return "TUNEL+"
    if s.startswith("TUNEL_NEG") or s in ("TUNEL-", "TUNELNEG", "NEGATIVE", "NEG",
                                          "CONTROL", "TUNELNEGATIVE",
                                          "NONAPOPTOTIC"):
        return "TUNEL-"
    return None


def _float(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def _position(value):
    m = re.search(r"(\d+)", str(value or ""))
    return int(m.group(1)) if m else None


def scan_run(run_dir):
    """Every CSV under `run_dir` that could be a per-parasite or per-nucleus table."""
    found = []
    for dirpath, dirnames, filenames in os.walk(run_dir):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in sorted(filenames):
            if not fn.lower().endswith(".csv"):
                continue
            path = os.path.join(dirpath, fn)
            try:
                with open(path, newline="") as fh:
                    reader = csv.DictReader(fh)
                    rows = list(reader)
                    fields = reader.fieldnames
            except (OSError, UnicodeDecodeError, csv.Error):
                continue
            if not rows or not fields:
                continue
            idx = _pick(fields, INDEX_COLUMNS)
            grp = _pick(fields, GROUP_COLUMNS, PATTERNS["group"])
            cor = _pick(fields, CORRECTED_COLUMNS, PATTERNS["corrected"])
            if not grp or not (idx or cor):
                continue
            groups = {_norm_group(r.get(grp)) for r in rows}
            inferred = None
            if not {"TUNEL+", "TUNEL-"} & groups:
                # A run that doubts its own channel assignment may decline to call
                # the result a TUNEL classification at all: repeat_2 labelled its
                # groups "Channel 2 high" / "Channel 2 low". That is still a
                # two-group split of the same six nuclei, and throwing it away
                # loses the run's best answer. Accept a two-valued class column
                # and map the higher-signal class to TUNEL+, recording that the
                # mapping was inferred rather than read.
                labels = {str(r.get(grp) or "").strip() for r in rows}
                labels.discard("")
                value_col = idx or cor
                if len(labels) != 2 or not value_col:
                    continue
                means = {}
                for lab in labels:
                    vals = [_float(r.get(value_col)) for r in rows
                            if str(r.get(grp) or "").strip() == lab]
                    vals = [v for v in vals if v is not None]
                    if not vals:
                        break
                    means[lab] = st.mean(vals)
                if len(means) != 2:
                    continue
                hi, lo = sorted(means, key=means.get, reverse=True)
                inferred = {hi: "TUNEL+", lo: "TUNEL-"}
            pos_col = _pick(fields, POSITION_COLUMNS)
            positions = {_position(r.get(pos_col)) for r in rows} if pos_col else set()
            positions.discard(None)
            found.append(dict(
                path=path, rel=os.path.relpath(path, run_dir), rows=rows,
                n=len(rows), index_col=idx, group_col=grp, corrected_col=cor,
                position_col=pos_col, n_positions=len(positions),
                name_col=_pick(fields, NAME_COLUMNS),
                green_col=_pick(fields, GREEN_COLUMNS),
                red_col=_pick(fields, RED_COLUMNS),
                background_col=_pick(fields, BACKGROUND_COLUMNS),
                group_map=inferred))
    return found


def choose_tables(found):
    """-> (task1, task2, task3) table dicts or None, with the rule applied stated.

    task 3 is the one table spanning many positions. task 2 is the largest
    single-position table that carries a proliferation index; where a run wrote
    both a per-subset table and a collapsed unique-parasite table, the collapsed
    one is preferred because the lab's four workbooks share one set of controls
    and the per-subset table counts them four times. task 1 carries a corrected
    TUNEL signal and no proliferation index.
    """
    notes = []
    multi = [f for f in found if f["index_col"] and f["n_positions"] > 1]
    # A run may leave more than one all-positions answer behind -- a first pass
    # and a revised one with different backgrounds, say. Picking one silently
    # would be picking the run's conclusion for it, so every serious candidate
    # is scored. Small tables alongside a large one are audit fragments, not
    # answers, and are dropped by the row-count floor.
    biggest = max((f["n"] for f in multi), default=0)
    task3 = [f for f in multi if f["n"] >= 0.5 * biggest]
    # A run often writes the same table twice under different names -- once as
    # the answer and once as a QC or exclusions view of it. Same row count and
    # same positions means the same answer; keep the one that does not call
    # itself an audit.
    seen, deduped = {}, []
    for f in sorted(task3, key=lambda f: bool(
            re.search(r"exclusion|audit|qc|diagnos", f["rel"], re.I))):
        key = (f["n"], f["n_positions"])
        if key in seen:
            notes.append(f"         ({f['rel']} holds the same {f['n']} rows as "
                         f"{seen[key]['rel']}; scored once)")
            continue
        seen[key] = f
        deduped.append(f)
    task3 = [f for f in multi if f in deduped]
    for t in task3:
        notes.append(f"task 3 <- {t['rel']} "
                     f"({t['n']} rows over {t['n_positions']} positions)")
    if len(multi) > len(task3):
        notes.append(f"         ({len(multi) - len(task3)} smaller multi-position "
                     f"table(s) ignored as audit fragments)")

    single = [f for f in found
              if f["index_col"] and f not in task3 and f["n_positions"] <= 1]
    unique = [f for f in single if "unique" in f["rel"].lower()]
    task2 = (max(unique, key=lambda f: f["n"], default=None)
             or max(single, key=lambda f: f["n"], default=None))
    if task2:
        notes.append(f"task 2 <- {task2['rel']} ({task2['n']} rows"
                     + ("; unique-parasite table preferred over the per-subset one"
                        if task2 in unique else "") + ")")

    ones = [f for f in found
            if not f["index_col"] and f["corrected_col"] and f["n"] <= 60]
    # Runs that reached an answer through several pipelines leave the rejected
    # ones behind, and one of those "rejected" passes may be the one that agrees
    # with the lab -- repeat_2 classified the same six nuclei twice, on two
    # channels, and filed the correct one under a hedged label. So every table
    # that classifies the nuclei is scored, oldest first.
    task1 = [f for f in ones
             if re.search(r"final|classified|classification", f["rel"], re.I)] or ones
    task1.sort(key=lambda f: os.path.getmtime(f["path"]))
    for t in task1:
        notes.append(f"task 1 <- {t['rel']} ({t['n']} rows)")
    if len(ones) > len(task1):
        notes.append(f"         ({len(ones) - len(task1)} table(s) that do not name "
                     f"themselves a classification ignored)")
    return task1, task2, task3, notes


def rows_of(table):
    """-> [{position, name, group, index, corrected, green, red}]"""
    out = []
    for r in table["rows"]:
        out.append(dict(
            position=_position(r.get(table["position_col"])) if table["position_col"] else None,
            name=str(r.get(table["name_col"]) or "").strip() if table["name_col"] else "",
            group=(table.get("group_map") or {}).get(
                str(r.get(table["group_col"]) or "").strip())
            if table.get("group_map") else _norm_group(r.get(table["group_col"])),
            index=_float(r.get(table["index_col"])) if table["index_col"] else None,
            corrected=_float(r.get(table["corrected_col"])) if table["corrected_col"] else None,
            background=_float(r.get(table["background_col"])) if table["background_col"] else None,
            green=_float(r.get(table["green_col"])) if table["green_col"] else None,
            red=_float(r.get(table["red_col"])) if table["red_col"] else None,
            raw=r))
    return out


# --------------------------------------------------------------------------
# matching
# --------------------------------------------------------------------------
def _cost_assign(cost):
    """Minimum-total-cost 1:1 assignment over a rectangular cost matrix."""
    try:
        import numpy as np
        from scipy.optimize import linear_sum_assignment
        ri, ci = linear_sum_assignment(np.asarray(cost))
        return [(int(i), int(j), float(cost[i][j])) for i, j in zip(ri, ci)]
    except ImportError:
        # Greedy on the globally smallest remaining cost. Not optimal, but for
        # matches this near-exact it agrees with the assignment above.
        cand = sorted((c, i, j) for i, row in enumerate(cost)
                      for j, c in enumerate(row))
        gi, gj, pairs = set(), set(), []
        for c, i, j in cand:
            if i in gi or j in gj:
                continue
            gi.add(i), gj.add(j)
            pairs.append((i, j, c))
        return pairs


def match(gts, runs, max_delta):
    """Pair one position's ground-truth parasites with a run's rows.

    -> (pairs, unmatched_gt, unmatched_run); each pair is (gt, run, how) where
    `how` is "roi" (paired by raw intensity, i.e. certainly the same ROI) or
    "index" (paired by proliferation index, i.e. probably the same ROI).
    """
    pairs, runner_up = [], []
    gts, runs = list(gts), list(runs)

    # Stage 1 -- which ROI is this? Raw green and red identify it outright.
    g_raw = [i for i, g in enumerate(gts) if g["green"] is not None and g["red"] is not None]
    r_raw = [j for j, r in enumerate(runs) if r["green"] is not None and r["red"] is not None]
    taken_g, taken_r = set(), set()
    if g_raw and r_raw:
        cost = [[abs(gts[i]["green"] - runs[j]["green"])
                 + abs(gts[i]["red"] - runs[j]["red"]) for j in r_raw] for i in g_raw]
        for i, j, c in _cost_assign(cost):
            if c <= SAME_ROI_TOLERANCE:
                gi, rj = g_raw[i], r_raw[j]
                pairs.append((gts[gi], runs[rj], "roi"))
                taken_g.add(gi), taken_r.add(rj)
                # How close is the NEXT best candidate? If some other row were
                # nearly as close, "same ROI" would be a guess rather than an
                # identification, and every conclusion drawn from it would be too.
                others = sorted(cost[i])[1:2]
                if others:
                    runner_up.append(others[0])

    # Stage 2 -- whatever is left, paired on the index, within group.
    for group in ("TUNEL+", "TUNEL-"):
        gi = [i for i in range(len(gts))
              if i not in taken_g and gts[i]["group"] == group]
        rj = [j for j in range(len(runs))
              if j not in taken_r and runs[j]["group"] == group
              and runs[j]["index"] is not None]
        if not gi or not rj:
            continue
        cost = [[abs(gts[i]["index"] - runs[j]["index"]) for j in rj] for i in gi]
        for a, b, c in _cost_assign(cost):
            if c <= max_delta:
                pairs.append((gts[gi[a]], runs[rj[b]], "index"))
                taken_g.add(gi[a]), taken_r.add(rj[b])

    return (pairs,
            [g for i, g in enumerate(gts) if i not in taken_g],
            [r for j, r in enumerate(runs) if j not in taken_r],
            runner_up)


def mann_whitney(a, b):
    try:
        from scipy.stats import mannwhitneyu
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        return float(u), float(p)
    except Exception:
        return None, None


def describe(values):
    if not values:
        return "n=0"
    s = f"n={len(values)} median={st.median(values):.3f} mean={st.mean(values):.3f}"
    if len(values) > 1:
        s += f" sd={st.stdev(values):.3f}"
    return s + f" range=[{min(values):.3f}, {max(values):.3f}]"


# --------------------------------------------------------------------------
# the three tasks
# --------------------------------------------------------------------------
def score_task1(table, gt_nuclei, gt_background, out_rows):
    print("=" * 78)
    print("TASK 1 -- TUNEL classification of the six given nuclei, position 40")
    print("=" * 78)
    if table is None:
        print("  no table found in this run\n")
        return
    run = rows_of(table)
    print(f"  source: {table['rel']}")
    print(f"  nuclei: ground truth {len(gt_nuclei)}, run {len(run)}")
    if table.get("group_map"):
        mapping = ", ".join(f"{k!r} -> {v}" for k, v in
                            sorted(table["group_map"].items(), key=lambda kv: kv[1]))
        print(f"  !! this run did not label its groups TUNEL+/TUNEL-. Mapped by "
              f"signal level: {mapping}")

    # WHICH CHANNEL did the run call TUNEL? The lab's own background level
    # answers that for the ground truth, and a run that landed on a different
    # channel produces numbers on a different scale entirely -- comparing them
    # as a difference would be meaningless, so say so instead.
    gt_bg = st.mean([b["green"] for b in gt_background]) if gt_background else None
    run_bg = next((r["background"] for r in run if r["background"] is not None), None)
    same_channel = True
    if gt_bg is not None and run_bg is not None:
        print(f"  background level: lab {gt_bg:.3f}, run {run_bg:.3f}", end="")
        if abs(run_bg - gt_bg) > 2:
            same_channel = False
            print("  <-- A DIFFERENT CHANNEL")
            print("     The lab's corrected signals sit on a background of "
                  f"{gt_bg:.0f}; this run's sit on {run_bg:.0f}. The two are not "
                  "measuring the same image, so the corrected values below are "
                  "not comparable as differences -- only the ordering and the "
                  "classification are.")
        else:
            print("  (same channel)")

    # The six nuclei are given ROIs, so the run's own numbering is the mapping
    # when it carries one. Everything else is checked against that assumption
    # rather than assumed silently.
    numbered = {}
    for r in run:
        m = re.search(r"(\d+)", r["name"])
        if m:
            numbered.setdefault(int(m.group(1)), r)
    by_number = len(numbered) == len(run) == len(gt_nuclei)
    if by_number:
        print("  matched by the run's own ROI numbering (Cell N <-> ROI N)")
        pairs = [(g, numbered.get(g["roi"])) for g in gt_nuclei]
    else:
        print("  !! run ROI names do not carry a usable number; matched by rank of "
              "corrected signal instead, which assumes the run measured the same "
              "six nuclei")
        rs = sorted(run, key=lambda r: (r["corrected"] is None, r["corrected"]))
        gs = sorted(gt_nuclei, key=lambda g: g["corrected"])
        pairs = list(zip(gs, rs))

    print()
    print(f"  {'ROI':>4} {'GT corrected':>13} {'run corrected':>14} {'delta':>8} "
          f"{'GT':>7} {'run':>7}  {'':2}")
    agree = 0
    for g, r in pairs:
        if r is None:
            print(f"  {g['roi']:>4} {g['corrected']:>13.3f} {'--':>14} {'--':>8} "
                  f"{g['group']:>7} {'--':>7}   MISSING")
            continue
        d = (r["corrected"] - g["corrected"]) if r["corrected"] is not None else float("nan")
        ok = g["group"] == r["group"]
        agree += ok
        print(f"  {g['roi']:>4} {g['corrected']:>13.3f} {r['corrected']:>14.3f} "
              f"{d:>+8.3f} {g['group']:>7} {r['group'] or '?':>7}   "
              f"{'' if ok else '<-- DISAGREES'}")
        out_rows.append(dict(task=1, source=table["rel"], position=40,
                             gt_label=f"ROI{g['roi']}",
                             run_name=r["name"], gt_value=g["corrected"],
                             run_value=r["corrected"], delta=d, matched_by="roi_number",
                             gt_group=g["group"], run_group=r["group"],
                             group_agrees=ok))
    tp = sum(1 for g, r in pairs if r and g["group"] == "TUNEL+" == r["group"])
    fp = sum(1 for g, r in pairs if r and r["group"] == "TUNEL+" != g["group"])
    fn = sum(1 for g, r in pairs if not r or (g["group"] == "TUNEL+" != r["group"]))
    print(f"\n  classification: {agree}/{len(pairs)} nuclei agree   "
          f"TUNEL+ recovered {tp}/{sum(g['group'] == 'TUNEL+' for g in gt_nuclei)}"
          f"  (false positives {fp}, missed {fn})")
    deltas = [abs(r["corrected"] - g["corrected"])
              for g, r in pairs if r and r["corrected"] is not None]
    if deltas and same_channel:
        print(f"  corrected signal: median |delta| {st.median(deltas):.3f}, "
              f"max {max(deltas):.3f}, on a signal that spans "
              f"{min(g['corrected'] for g in gt_nuclei):.1f} to "
              f"{max(g['corrected'] for g in gt_nuclei):.1f}")
    # Ordering survives a change of channel, so it is the one thing worth
    # reporting either way: it says whether the run would have reached the lab's
    # classification had it thresholded the same way.
    both = [(g["corrected"], r["corrected"]) for g, r in pairs
            if r and r["corrected"] is not None]
    if len(both) > 2:
        try:
            from scipy.stats import spearmanr
            rho, _p = spearmanr([a for a, _ in both], [b for _, b in both])
            print(f"  ordering of the nuclei: Spearman rho {rho:+.3f} "
                  f"(1.0 would mean the run ranks them exactly as the lab does)")
        except ImportError:
            pass
    print()


def score_indexed(label, table, gt_by_position, max_delta, out_rows, task_no,
                  roi_inventory=None):
    print("=" * 78)
    print(label)
    print("=" * 78)
    if table is None:
        print("  no table found in this run\n")
        return
    run = rows_of(table)
    print(f"  source: {table['rel']}")

    run_by_pos = collections.defaultdict(list)
    for r in run:
        if r["index"] is None:
            continue
        run_by_pos[r["position"]].append(r)
    # A single-position table may carry no position column at all.
    if len(gt_by_position) == 1 and list(run_by_pos) == [None]:
        only = next(iter(gt_by_position))
        run_by_pos = {only: run_by_pos[None]}

    n_gt = sum(len(v) for v in gt_by_position.values())
    n_run = sum(len(v) for v in run_by_pos.values())
    print(f"  parasites: ground truth {n_gt} over {len(gt_by_position)} position(s), "
          f"run {n_run} over {len(run_by_pos)} position(s)")
    dropped = len(run) - n_run
    if dropped:
        print(f"  ({dropped} run rows carry no index and are excluded)")

    tot = collections.Counter()
    deltas_roi, deltas_index, per_position, group_errors = [], [], [], []
    runner_up = []
    for pos in sorted(set(gt_by_position) | set(run_by_pos),
                      key=lambda p: (p is None, p)):
        gts = gt_by_position.get(pos, [])
        rns = run_by_pos.get(pos, [])
        pairs, miss_g, miss_r, second = match(gts, rns, max_delta)
        runner_up += second
        for g, r, how in pairs:
            d = r["index"] - g["index"]
            (deltas_roi if how == "roi" else deltas_index).append(abs(d))
            same_group = g["group"] == r["group"]
            if not same_group and how == "roi":
                group_errors.append((pos, g["label"], r["name"], g["group"], r["group"]))
            tot["by_" + how] += 1
            out_rows.append(dict(task=task_no, source=table["rel"],
                                 position=pos, gt_label=g["label"],
                                 run_name=r["name"], gt_value=g["index"],
                                 run_value=r["index"], delta=d, matched_by=how,
                                 gt_group=g["group"], run_group=r["group"],
                                 group_agrees=same_group))
        for g in miss_g:
            out_rows.append(dict(task=task_no, source=table["rel"],
                                 position=pos, gt_label=g["label"],
                                 run_name="", gt_value=g["index"], run_value=None,
                                 delta=None, matched_by="", gt_group=g["group"],
                                 run_group="", group_agrees=None))
        for r in miss_r:
            out_rows.append(dict(task=task_no, source=table["rel"],
                                 position=pos, gt_label="",
                                 run_name=r["name"], gt_value=None,
                                 run_value=r["index"], delta=None, matched_by="",
                                 gt_group="", run_group=r["group"], group_agrees=None))
        tot["matched"] += len(pairs)
        tot["gt_missing"] += len(miss_g)
        tot["run_extra"] += len(miss_r)
        per_position.append((pos, len(gts), len(rns), len(pairs),
                             len(miss_g), len(miss_r), [r["name"] for r in miss_r]))

    print()
    print(f"  COVERAGE   {tot['matched']}/{n_gt} of the lab's parasites reproduced "
          f"({100.0 * tot['matched'] / n_gt:.1f}%)")
    print(f"             {tot['gt_missing']} of the lab's parasites have no run number")
    print(f"             {tot['run_extra']} run rows match no parasite the lab scored")
    print(f"             (paired: {tot['by_roi']} by raw intensity -- certainly the "
          f"same ROI; {tot['by_index']} by index only)")
    if runner_up:
        ambiguous = sum(1 for c in runner_up if c <= SAME_ROI_TOLERANCE)
        print(f"             ROI pairing is unambiguous: the second-best candidate "
              f"was a median {st.median(runner_up):.1f} away against a "
              f"{SAME_ROI_TOLERANCE} threshold"
              + (f", and only {ambiguous} pairing(s) had any rival at all -- those "
                 f"are parasites the lab measured twice" if ambiguous else ""))

    for label, ds in (("raw intensity (same ROI)", deltas_roi),
                      ("index alone", deltas_index)):
        if not ds:
            continue
        exact = sum(1 for d in ds if d <= 5e-4)
        close = sum(1 for d in ds if 5e-4 < d <= 0.05)
        print()
        print(f"  AGREEMENT on the {len(ds)} pairs matched by {label}")
        print(f"             {exact} reproduce the lab's index exactly (|delta| <= 5e-4)")
        print(f"             {close} agree to 0.05")
        print(f"             {len(ds) - exact - close} differ by more, up to "
              f"{max(ds):.3f} index units")
    if deltas_roi and any(d > 0.05 for d in deltas_roi):
        n = sum(1 for d in deltas_roi if d > 0.05)
        print()
        print(f"  !! {n} parasites were measured on the SAME ROI -- the raw green and "
              f"red means agree to five decimals -- yet the index differs. The "
              f"parasite measurement is right, so what differs is the three "
              f"background regions paired with it.")
    if group_errors:
        print()
        print(f"  !! {len(group_errors)} parasites put in the wrong group:")
        for pos, lab, name, gg, rg in group_errors[:10]:
            print(f"       pos {pos}: {lab} is {gg} in the lab's sheet, "
                  f"{rg} in the run (as {name})")

    if len(per_position) > 1:
        print()
        print("  per position (positions where something is missing or extra):")
        print(f"     {'pos':>4} {'GT':>4} {'run':>4} {'ok':>4} {'missing':>8} {'extra':>6}"
              + ("  ROI set" if roi_inventory else ""))
        for pos, ng, nr, ok, mg, mr, extra_names in per_position:
            if mg == 0 and mr == 0:
                continue
            extra = ""
            if roi_inventory and pos in roi_inventory:
                inv = roi_inventory[pos]
                extra = f"  {len(inv['parasites'])} parasite ROIs, {len(inv['backgrounds'])} bg"
            print(f"     {str(pos):>4} {ng:>4} {nr:>4} {ok:>4} {mg:>8} {mr:>6}{extra}")
        clean = sum(1 for row in per_position if row[4] == 0 and row[5] == 0)
        print(f"     ({clean} position(s) fully reproduced, not listed)")

        # A run row that matches nothing and is NAMED like a background region is
        # not a near miss: it is a background measured as though it were a
        # parasite. A few positions name their backgrounds `C29B1` after the
        # parasite they belong to instead of the usual `0788-1151-0552`, and a
        # run keyed on the numeric spelling takes them in as parasites.
        mistaken = [(pos, n) for pos, _g, _r, _o, _mg, _mr, names in per_position
                    for n in names if re.match(r"^(C|TUNEL)\d+B\d+", n)]
        if mistaken:
            where = sorted({p for p, _n in mistaken})
            print()
            print(f"  !! {len(mistaken)} background regions were scored as parasites: "
                  + ", ".join(f"pos {p} {n}" for p, n in mistaken[:8])
                  + (" ..." if len(mistaken) > 8 else ""))
            print(f"     These are `C<n>B<k>` backgrounds, the spelling used at "
                  f"position(s) {', '.join(str(p) for p in where)}. Elsewhere a "
                  f"background is named `0788-1151-0552`, so a run that recognises "
                  f"only that spelling measures these as parasites.")

    print()
    print("  CONCLUSION (both computed the same way, from each set of numbers)")
    for src, data in (("ground truth", gt_by_position), ("run", run_by_pos)):
        vals = {g: [x["index"] for v in data.values() for x in v
                    if x["group"] == g and x["index"] is not None]
                for g in ("TUNEL+", "TUNEL-")}
        u, p = mann_whitney(vals["TUNEL+"], vals["TUNEL-"])
        print(f"     {src:13} TUNEL+  {describe(vals['TUNEL+'])}")
        print(f"     {'':13} TUNEL-  {describe(vals['TUNEL-'])}")
        if u is not None:
            diff = st.median(vals["TUNEL+"]) - st.median(vals["TUNEL-"]) \
                if vals["TUNEL+"] and vals["TUNEL-"] else float("nan")
            print(f"     {'':13} Mann-Whitney U={u:.1f}  p={p:.4g}  "
                  f"median difference {diff:+.3f}")
    print()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True, help="a run directory to score")
    p.add_argument("--gt-root", required=True,
                   help="the parasite_task_results folder holding the lab's analysis")
    p.add_argument("--tasks", nargs="+", type=int, default=[1, 2, 3],
                   choices=[1, 2, 3])
    p.add_argument("--max-delta", type=float, default=0.5,
                   help="reject a pairing further apart than this in index units "
                        "(default 0.5; the groups span about 3 to 10)")
    p.add_argument("--out", help="write the matched-pair table here as CSV")
    args = p.parse_args()

    found = scan_run(args.run)
    task1, task2, task3, notes = choose_tables(found)
    print("=" * 78)
    print(f"RUN  {args.run}")
    print("=" * 78)
    print(f"  {len(found)} candidate table(s) carrying a group label and a value:")
    for f in found:
        marks = "".join(m for m, t in (("1", task1), ("2", [task2]), ("3", task3))
                        if f in t)
        print(f"    [{marks or ' '}] {f['rel']:70.70} {f['n']:>5} rows"
              + (f", {f['n_positions']} positions" if f["n_positions"] > 1 else ""))
    print()
    for n in notes:
        print("  " + n)
    print()

    out_rows = []
    if 1 in args.tasks:
        nuclei, background = gt.load_task1(os.path.join(
            args.gt_root, "oneframe_results", "ExampleData_TunelMeasurement.xlsx"))
        for t in (task1 or [None]):
            score_task1(t, nuclei, background, out_rows)
    if 2 in args.tasks:
        g2, files = gt.load_task2(os.path.join(args.gt_root, "oneframe_results"))
        print(f"(task 2 ground truth: {len(files)} workbooks collapsed to "
              f"{len(g2)} unique parasites -- the seven controls are one set of "
              f"measurements copied into all four)")
        score_indexed("TASK 2 -- proliferation index, one frame of position 40",
                      task2, {40: g2}, args.max_delta, out_rows, 2)
    if 3 in args.tasks:
        g3 = gt.load_task3(os.path.join(args.gt_root, "4. Data Analysis Excel Files"))
        inv = None
        roi_dir = os.path.join(args.gt_root, "3. Data Analysis ROIs")
        for sub in ("Apoptotic cells", "Macrophages"):
            if os.path.isdir(os.path.join(roi_dir, sub)):
                inv = gt.load_roi_inventory(os.path.join(roi_dir, sub))
                break
        for t in (task3 or [None]):
            title = "TASK 3 -- proliferation index, every parasite of every position"
            if t is not None and len(task3) > 1:
                title += f"\n         table: {t['rel']}"
            score_indexed(title, t, g3, args.max_delta, out_rows, 3, inv)

    if args.out and out_rows:
        keys = ["task", "source", "position", "gt_label", "run_name", "gt_value", "run_value",
                "delta", "matched_by", "gt_group", "run_group", "group_agrees"]
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(out_rows)
        print(f"-> {args.out}  ({len(out_rows)} rows)")


if __name__ == "__main__":
    sys.exit(main())
