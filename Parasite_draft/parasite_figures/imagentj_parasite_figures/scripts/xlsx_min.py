#!/usr/bin/env python3
"""Minimal .xlsx reader — no openpyxl needed.

An xlsx is a zip of XML, and everything these ground-truth workbooks use is
shallow: shared strings, inline strings, numbers, and cached formula results.
The "10-R/G" column is a formula, so what is read back is the value Excel
cached the last time it recalculated; a workbook saved by a tool that does not
cache would give None there, which is reported rather than silently treated as
a missing parasite.

    read(path) -> {sheet_name: [[cell, ...], ...]}   cells are str | float | None
"""

import re
import zipfile
import xml.etree.ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RELNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _col_index(ref):
    """'C7' -> 2. Cells with no value are skipped in the XML, so the column
    letter is the only thing that says where a cell actually sits."""
    m = re.match(r"([A-Z]+)", ref or "")
    if not m:
        return 0
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read(path):
    z = zipfile.ZipFile(path)
    names = set(z.namelist())

    shared = []
    if "xl/sharedStrings.xml" in names:
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(f"{NS}si"):
            shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))

    rels = {}
    if "xl/_rels/workbook.xml.rels" in names:
        for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels")):
            rels[r.get("Id")] = r.get("Target")

    sheets = []
    for sh in ET.fromstring(z.read("xl/workbook.xml")).iter(f"{NS}sheet"):
        tgt = rels.get(sh.get(f"{RELNS}id"), "")
        if tgt and not tgt.startswith("xl/"):
            tgt = "xl/" + tgt.lstrip("/")
        sheets.append((sh.get("name"), tgt))

    out = {}
    for name, target in sheets:
        if target not in names:
            continue
        rows = []
        for row in ET.fromstring(z.read(target)).iter(f"{NS}row"):
            cells = []
            for c in row.findall(f"{NS}c"):
                i = _col_index(c.get("r"))
                while len(cells) < i:
                    cells.append(None)
                t = c.get("t")
                if t == "inlineStr":
                    v = "".join(x.text or "" for x in c.iter(f"{NS}t"))
                else:
                    ve = c.find(f"{NS}v")
                    v = ve.text if ve is not None else None
                    if v is not None:
                        if t == "s":
                            v = shared[int(v)]
                        elif t in (None, "n"):
                            try:
                                v = float(v)
                            except ValueError:
                                pass
                cells.append(v)
            rows.append(cells)
        out[name] = rows
    return out


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print("=" * 76)
        print(p)
        for sheet, rows in read(p).items():
            print(f"-- sheet {sheet!r}: {len(rows)} rows")
            for i, r in enumerate(rows):
                cells = ["" if c is None else
                         (f"{c:.6g}" if isinstance(c, float) else str(c)) for c in r]
                if any(cells):
                    print(f"{i:4d} | " + " | ".join(f"{c:>12}" for c in cells))
