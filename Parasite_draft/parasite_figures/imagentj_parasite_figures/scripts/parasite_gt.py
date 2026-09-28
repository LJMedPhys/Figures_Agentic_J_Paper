#!/usr/bin/env python3
"""Ground truth for the three Leishmania / parasite tasks, as the lab recorded it.

Three sources, two layouts:

  task 1  oneframe_results/ExampleData_TunelMeasurement.xlsx
          Six nuclei ROIs of position 40 with their raw TUNEL mean, the shared
          background ROIs below them, a `Corrected green:` column, and a
          `TUNEL+` flag column carrying a 1 on the cells the lab called positive.

  task 2  oneframe_results/EvaluationCalculation{1..4}.xlsx
          One workbook per parasite subset (set_1..set_4). Each holds TWO
          apoptotic parasites, unique to that subset, and the SAME seven control
          parasites, copied into all four workbooks. Pooling the four files
          naively gives 36 parasites where there are only 15 -- eight TUNEL+ and
          seven TUNEL-, the controls being one set of measurements repeated.
          `unique=True` collapses them and is the default.

  task 3  "4. Data Analysis Excel Files"/Position<N>_<k>.xlsx
          The same summary layout, `k` files per position, nine parasites per
          file (two `Apoptotic`, seven `Control`). Here nothing repeats: every
          row across a position's files is a distinct parasite.

TWO TRAPS, both real and both silent if ignored
-----------------------------------------------
1. `Control<n>` in the summary block is a TUNEL-NEGATIVE HOST, i.e. a parasite
   living in a cell that did not die. `Control area <n>` in the ROI-order legend
   further down the same sheet is a BACKGROUND REGION. They are unrelated. The
   ROI files spell the first `C<n>.roi` and the second `0788-1151-0552.roi`.

2. The raw intensity columns are NOT on a consistent scale. These sheets were
   made in a German-locale Excel and some were pasted from ImageJ as text, so
   the decimal comma was lost: `281,880` survived as `281.880` in one file and
   became the integer `310010` in another. The ratio columns `R/G` and `10-R/G`
   were computed by Excel from those numbers and are scale-free, so THEY are
   correct everywhere. Every intensity this module returns therefore carries a
   `scale` field saying which power of ten was divided out, and matching against
   a run is done on the index, never on the raw intensities.
"""

import glob
import math
import os
import re

import xlsx_min

# Summary-block headers, spelled as the workbooks actually spell them (several
# are typos that survived copy-paste: "Red GB1", "RedBG3").
SUMMARY_COLUMNS = {
    "green": "green", "red": "red",
    "corrected green": "cgreen", "corrected red": "cred",
    "r/g": "rg", "10-r/g": "index",
}
# What a plausible 8/16-bit-ish ROI mean looks like, used only to undo the lost
# decimal comma. A parasite mean of 310010 is not an intensity; 310.010 is.
PLAUSIBLE_MEAN = (1.0, 70000.0)


def _norm(s):
    return re.sub(r"\s+", " ", str(s)).strip().lower()


def _descale(value):
    """Undo a lost decimal comma. Returns (value, power_of_ten_removed)."""
    if value is None:
        return None, 0
    k = 0
    while abs(value) > PLAUSIBLE_MEAN[1] and k < 6:
        value /= 1000.0
        k += 3
    return value, k


def read_summary_block(path):
    """Parse one workbook that uses the `Green | ... | 10-R/G` summary layout."""
    out = []
    for sheet, rows in xlsx_min.read(path).items():
        header_row = None
        for i, r in enumerate(rows):
            if any(c is not None and _norm(c) == "10-r/g" for c in r):
                header_row = i
                break
        if header_row is None:
            continue
        col, bg_col = {}, {"green": [], "red": []}
        for j, c in enumerate(rows[header_row]):
            n = _norm(c) if c is not None else ""
            key = SUMMARY_COLUMNS.get(n)
            if key and key not in col:
                col[key] = j
            # The three background columns per channel, spelled inconsistently:
            # "Green BG1", "Red GB1", "RedBG3" all occur.
            m = re.match(r"^(green|red)\s*(bg|gb)\s*\d$", n)
            if m:
                bg_col[m.group(1)].append(j)
        if "index" not in col or "green" not in col:
            continue
        # The label sits in the column immediately left of `Green`; that header
        # cell is blank, so it cannot be found by name.
        label_col = col["green"] - 1
        for r in rows[header_row + 1:]:
            if len(r) <= label_col:
                continue
            label = r[label_col]
            if not isinstance(label, str) or not label.strip():
                continue
            m = re.match(r"^\s*(apoptotic|nonapoptotic|control)\s*(\d*)\s*$", label, re.I)
            if not m:
                continue

            def cell(key):
                j = col.get(key)
                v = r[j] if j is not None and j < len(r) else None
                return v if isinstance(v, (int, float)) else None

            # The lost decimal comma is a property of the ROW, not of each cell:
            # every number pasted from the same ImageJ table lost it together. A
            # corrected value is one to two orders of magnitude smaller than the
            # raw mean it came from, so descaling it on its own magnitude leaves
            # it untouched while the raw means are fixed -- and the two columns
            # then disagree by a factor of a thousand. So the factor is taken
            # from the raw means, where it is detectable, and applied to the rest.
            green, sg = _descale(cell("green"))
            red, sr = _descale(cell("red"))
            factor = 10.0 ** max(sg, sr)
            cgreen = cell("cgreen") / factor if cell("cgreen") is not None else None
            cred = cell("cred") / factor if cell("cred") is not None else None
            index = cell("index")
            derived = False
            if index is None and cred is not None and cgreen:
                index, derived = 10 - cred / cgreen, True
            if index is None:
                continue
            def bg(channel):
                vals = []
                for j in bg_col[channel]:
                    v = r[j] if j < len(r) else None
                    if isinstance(v, (int, float)):
                        vals.append(v / factor)
                return vals

            out.append(dict(
                source=os.path.basename(path), sheet=sheet, label=label.strip(),
                number=int(m.group(2)) if m.group(2) else None,
                group="TUNEL+" if m.group(1).lower() == "apoptotic" else "TUNEL-",
                green=green, red=red, cgreen=cgreen, cred=cred,
                bg_green=bg("green"), bg_red=bg("red"),
                index=index, index_derived=derived, scale=max(sg, sr)))
    return out


# --------------------------------------------------------------------------
# task 1 -- TUNEL classification of the six given nuclei of position 40
# --------------------------------------------------------------------------
def load_task1(path):
    """-> (nuclei, background_rois). Nuclei carry raw/corrected green and group."""
    rows = next(iter(xlsx_min.read(path).values()))
    header = None
    for i, r in enumerate(rows):
        if any(c is not None and _norm(c) == "roi" for c in r):
            header = i
            break
    if header is None:
        raise ValueError(f"no `ROI` header in {path}")
    # `TUNEL+` appears twice in this header: once as the flag column carrying a
    # 1, and again further right over the split-out positive/negative values.
    # First occurrence wins, which is the flag.
    fields = {"roi": "roi", "area": "area", "green": "green", "x": "x", "y": "y",
              "corrected green:": "corrected", "tunel+": "flag"}
    col = {}
    for j, c in enumerate(rows[header]):
        key = fields.get(_norm(c) if c is not None else "")
        if key and key not in col:
            col[key] = j
    nuclei, background, in_background = [], [], False
    for r in rows[header + 1:]:
        first = next((c for c in r if c is not None and str(c).strip() != ""), None)
        if isinstance(first, str) and _norm(first) == "background":
            in_background = True
            continue

        def cell(key):
            j = col.get(key)
            v = r[j] if j is not None and j < len(r) else None
            return v if isinstance(v, (int, float)) else None

        if cell("roi") is None:
            continue
        rec = dict(roi=int(cell("roi")), area=cell("area"), green=cell("green"),
                   x=cell("x"), y=cell("y"))
        if in_background:
            background.append(rec)
        else:
            rec["corrected"] = cell("corrected")
            rec["group"] = "TUNEL+" if cell("flag") else "TUNEL-"
            nuclei.append(rec)
    return nuclei, background


# --------------------------------------------------------------------------
# task 2 -- proliferation index of the parasites in one frame of position 40
# --------------------------------------------------------------------------
def load_task2(folder, unique=True):
    files = sorted(glob.glob(os.path.join(folder, "EvaluationCalculation*.xlsx")),
                   key=lambda p: int(re.search(r"(\d+)\.xlsx$", p).group(1)))
    files = [f for f in files if not os.path.basename(f).startswith((".", "~$"))]
    rows = []
    for k, p in enumerate(files, start=1):
        for r in read_summary_block(p):
            r["subset"] = f"set_{k}"
            rows.append(r)
    if not unique:
        return rows, files
    # The seven controls are one measurement pasted into all four workbooks:
    # identical label AND identical index. Collapse those, keep everything else.
    seen, out = {}, []
    for r in rows:
        key = (r["group"], r["label"], round(r["index"], 9))
        if key in seen:
            seen[key]["subsets"].append(r["subset"])
            continue
        r = dict(r, subsets=[r["subset"]])
        seen[key] = r
        out.append(r)
    return out, files


# --------------------------------------------------------------------------
# task 3 -- every parasite of every position
# --------------------------------------------------------------------------
def load_task3(folder):
    """-> {position:int -> [parasite, ...]}, each carrying source workbook."""
    by_position = {}
    for p in sorted(glob.glob(os.path.join(folder, "Position*.xlsx"))):
        base = os.path.basename(p)
        if base.startswith((".", "~$")):
            continue
        m = re.match(r"Position(\d+)_(\d+)\.xlsx$", base)
        if not m:
            raise ValueError(f"unexpected workbook name: {base}")
        pos, part = int(m.group(1)), int(m.group(2))
        for r in read_summary_block(p):
            r["position"], r["part"] = pos, part
            by_position.setdefault(pos, []).append(r)
    return by_position


# --------------------------------------------------------------------------
# the ROI sets, which say which parasites EXIST even where no workbook scored
# them, and which name the background regions
# --------------------------------------------------------------------------
def load_roi_inventory(folder):
    """-> {position:int -> {'parasites': [...], 'backgrounds': [...]}}.

    `C<n>` and `TUNEL<n>` are parasites. Backgrounds come in two spellings: the
    ImageJ default `0788-1151-0552` in most positions, and `C29B1` in positions
    39 and 40, where whoever drew them named the background after its parasite.
    A run that assumes only the first spelling scores those as parasites.
    """
    out = {}
    for d in sorted(os.listdir(folder)):
        m = re.match(r"Position_(\d+)$", d)
        if not m or not os.path.isdir(os.path.join(folder, d)):
            continue
        parasites, backgrounds = [], []
        for f in sorted(os.listdir(os.path.join(folder, d))):
            if not f.endswith(".roi"):
                continue
            name = f[:-4]
            if re.match(r"^(C|TUNEL)\d+B\d+(-\d+)?$", name) or name[0].isdigit():
                backgrounds.append(name)
            elif re.match(r"^(C|TUNEL)\d+", name):
                parasites.append(name)
        out[int(m.group(1))] = dict(parasites=parasites, backgrounds=backgrounds)
    return out


# --------------------------------------------------------------------------
# WHERE the lab's objects are, for runs that segmented their own
# --------------------------------------------------------------------------
# A human-in-the-loop run measured the lab's ROIs, so its raw intensities
# identify which parasite each row is. An auto-mode run segmented its own
# objects, so nothing in its output matches the lab's numbers and the only way
# to ask "did it find this parasite?" is position. The workbooks hold no
# coordinates; the .roi files do.
PIXEL_UM = 0.32231030776746455


def load_task1_coords(one_frame_dir, nuclei=None):
    """The six nuclei of position 40, in FULL-image pixels.

    `nuclei_tunel_roi/Cell N.roi` is stored relative to `TUNEL_Overview.roi`, so
    the overview origin has to be added. Getting this wrong by the overview
    offset is the mistake one HIL run made and had to be told about, so when the
    workbook rows are passed in, their calibrated X/Y is used to CHECK the
    result rather than trusting the convention.
    """
    import ij_roi
    origin = ij_roi.read(os.path.join(one_frame_dir, "TUNEL_Overview.roi"))
    rois = ij_roi.read_dir(os.path.join(one_frame_dir, "nuclei_tunel_roi"))
    out = []
    for r in rois:
        m = re.match(r"^Cell\s*(\d+)$", r["name"])
        rec = dict(name=r["name"], roi=int(m.group(1)) if m else None,
                   kind="nucleus" if m else "background",
                   x=r["cx"] + origin["left"], y=r["cy"] + origin["top"],
                   w=r["w"], h=r["h"])
        out.append(rec)
    nucleus_rois = [r for r in out if r["kind"] == "nucleus"]
    check = None
    if nuclei:
        # The workbook's X/Y are micrometres inside the crop; the ROI centre is
        # pixels inside the crop. If the two disagree the origin is wrong.
        offs = []
        by_number = {n["roi"]: n for n in nuclei}
        for r in nucleus_rois:
            n = by_number.get(r["roi"])
            if not n or n.get("x") is None:
                continue
            offs.append(math.hypot((n["x"] / PIXEL_UM + origin["left"]) - r["x"],
                                   (n["y"] / PIXEL_UM + origin["top"]) - r["y"]))
        if offs:
            check = max(offs)
    return nucleus_rois, [r for r in out if r["kind"] == "background"], check


def load_task2_coords(one_frame_dir):
    """The unique parasites of the position-40 frame, in FULL-image pixels.

    `parasite_roi/set_N/` holds `C1..C7` (the same seven controls in every set,
    identical geometry) plus that set's two TUNEL-positive parasites. The
    controls are deduplicated on position; `C<n>` corresponds to the workbooks'
    `Control<n>`, which is how a matched detection gets a ground-truth index.
    """
    import ij_roi
    rois = ij_roi.read_dir(os.path.join(one_frame_dir, "parasite_roi"))
    out, seen = [], {}
    for r in rois:
        name = r["name"]
        if re.match(r"^\d{4}-\d{4}-\d{4}$", name):
            continue                      # a background region
        control = re.match(r"^C(\d+)$", name)
        key = (name, round(r["cx"], 1), round(r["cy"], 1))
        if control and key in seen:
            continue                      # the same control, copied into another set
        seen[key] = True
        out.append(dict(name=name, subset=r.get("folder", ""),
                        group="TUNEL-" if control else "TUNEL+",
                        gt_label=f"Control{control.group(1)}" if control else None,
                        x=r["cx"], y=r["cy"], w=r["w"], h=r["h"]))
    return out


if __name__ == "__main__":
    import sys
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    nuclei, bg = load_task1(os.path.join(
        root, "oneframe_results", "ExampleData_TunelMeasurement.xlsx"))
    print(f"task 1: {len(nuclei)} nuclei "
          f"({sum(n['group'] == 'TUNEL+' for n in nuclei)} TUNEL+), "
          f"{len(bg)} background ROIs")
    t2, files = load_task2(os.path.join(root, "oneframe_results"))
    print(f"task 2: {len(t2)} unique parasites from {len(files)} workbooks "
          f"({sum(r['group'] == 'TUNEL+' for r in t2)} TUNEL+)")
    t3 = load_task3(os.path.join(root, "4. Data Analysis Excel Files"))
    n = sum(len(v) for v in t3.values())
    print(f"task 3: {n} parasites over {len(t3)} positions "
          f"({sum(r['group'] == 'TUNEL+' for v in t3.values() for r in v)} TUNEL+)")
    inv = load_roi_inventory(os.path.join(root, "3. Data Analysis ROIs", "Apoptotic cells"))
    print(f"ROI sets: {sum(len(v['parasites']) for v in inv.values())} parasite ROIs, "
          f"{sum(len(v['backgrounds']) for v in inv.values())} background ROIs, "
          f"over {len(inv)} positions")
