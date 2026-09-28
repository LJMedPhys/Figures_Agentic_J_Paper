#!/usr/bin/env python3
"""Minimal multi-sheet .xlsx writer — no openpyxl, matching xlsx_min.py's reader.

    write("out.xlsx", [
        ("Sheet name", [["header", "row"], ["data", 1.0], ...]),
        ...
    ])

Cells are str, int, float, bool or None. Strings are written inline, so there is
no shared-string table to keep consistent. The first row of each sheet is styled
bold and frozen, and column widths are set from the content, because a results
table nobody can read without resizing forty columns is not a deliverable.

A cell whose string starts with `=` is written as a formula, and a string that
looks like a number is still written as text, so an ID like `007` survives.
"""

import re
import zipfile

MAX_WIDTH = 60


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def _col_letter(i):
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _cell(ref, value, style):
    st = f' s="{style}"' if style else ""
    if value is None or value == "":
        return f'<c r="{ref}"{st}/>'
    if isinstance(value, bool):
        return f'<c r="{ref}"{st} t="b"><v>{1 if value else 0}</v></c>'
    if isinstance(value, (int, float)):
        return f'<c r="{ref}"{st}><v>{value!r}</v></c>'
    text = str(value)
    if text.startswith("="):
        return f'<c r="{ref}"{st}><f>{_esc(text[1:])}</f></c>'
    return (f'<c r="{ref}"{st} t="inlineStr"><is><t xml:space="preserve">'
            f'{_esc(text)}</t></is></c>')


def _sheet_xml(rows, header=True):
    widths = {}
    for row in rows:
        for j, v in enumerate(row):
            n = len(str(v)) if v is not None else 0
            widths[j] = min(MAX_WIDTH, max(widths.get(j, 8), n + 2))
    cols = "".join(f'<col min="{j + 1}" max="{j + 1}" width="{w}" customWidth="1"/>'
                   for j, w in sorted(widths.items()))
    body = []
    for i, row in enumerate(rows):
        style = 1 if (header and i == 0) else 0
        cells = "".join(_cell(f"{_col_letter(j)}{i + 1}", v, style)
                        for j, v in enumerate(row))
        body.append(f'<row r="{i + 1}">{cells}</row>')
    pane = ('<sheetViews><sheetView workbookViewId="0">'
            '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
            '</sheetView></sheetViews>') if header else ""
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f'{pane}<cols>{cols}</cols><sheetData>{"".join(body)}</sheetData></worksheet>')


STYLES = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
          '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
          '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
          '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
          '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
          '<borders count="1"><border/></borders>'
          '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0"/></cellStyleXfs>'
          '<cellXfs count="2"><xf numFmtId="0" fontId="0" xfId="0"/>'
          '<xf numFmtId="0" fontId="1" xfId="0" applyFont="1"/></cellXfs>'
          '</styleSheet>')


def _safe_name(name, used):
    # Excel forbids : \ / ? * [ ] and caps sheet names at 31 characters.
    clean = re.sub(r"[:\\/?*\[\]]", "-", str(name))[:31] or "Sheet"
    base, k = clean, 2
    while clean in used:
        clean = f"{base[:28]}_{k}"
        k += 1
    used.add(clean)
    return clean


def write(path, sheets, header=True):
    """sheets: [(name, rows)] where rows is a list of lists."""
    used = set()
    named = [(_safe_name(n, used), r) for n, r in sheets]
    n = len(named)

    types = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
             '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
             '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
             '<Default Extension="xml" ContentType="application/xml"/>',
             '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
             '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>']
    for i in range(n):
        types.append(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" '
                     'ContentType="application/vnd.openxmlformats-officedocument.'
                     'spreadsheetml.worksheet+xml"/>')
    types.append("</Types>")

    wb_sheets = "".join(
        f'<sheet name="{_esc(name)}" sheetId="{i + 1}" r:id="rId{i + 1}"/>'
        for i, (name, _rows) in enumerate(named))
    workbook = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
                f'<sheets>{wb_sheets}</sheets></workbook>')
    wb_rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
               '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    for i in range(n):
        wb_rels.append(f'<Relationship Id="rId{i + 1}" Type="http://schemas.'
                       'openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                       f'Target="worksheets/sheet{i + 1}.xml"/>')
    wb_rels.append(f'<Relationship Id="rId{n + 1}" Type="http://schemas.'
                   'openxmlformats.org/officeDocument/2006/relationships/styles" '
                   'Target="styles.xml"/>')
    wb_rels.append("</Relationships>")

    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/'
                 'officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                 '</Relationships>')

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", "".join(types))
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", workbook)
        z.writestr("xl/_rels/workbook.xml.rels", "".join(wb_rels))
        z.writestr("xl/styles.xml", STYLES)
        for i, (_name, rows) in enumerate(named):
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml", _sheet_xml(rows, header))
    return path


if __name__ == "__main__":
    import xlsx_min
    p = write("/tmp/_xlsx_write_selftest.xlsx",
              [("First", [["name", "value", "flag"], ["alpha", 1.5, True],
                          ["beta", -2, False], ["gamma", None, None]]),
               ("Second: bad/name*", [["only", "one"], ["row", 42]])])
    back = xlsx_min.read(p)
    for sheet, rows in back.items():
        print(f"{sheet!r}: {rows}")
