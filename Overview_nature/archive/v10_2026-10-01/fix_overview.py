"""Make a Nature Methods-compliant copy of ../Overview.svg without changing its layout.

The drawing keeps its coordinates; only the physical size is set to 180 mm wide
(double column), so every length shrinks by K = 180 / 285.48. Font sizes are then
chosen so that they PRINT at the target point size after that shrink:

    size_in_user_units = pt * 0.3528 / K

What this script changes, and nothing else:
  * width/height -> 180 mm (viewBox untouched)
  * every text size, mapped by its role (see SIZE_MAP), Arial instead of generic 'sans'
  * empty leftover <text> elements removed
  * the hairline outline on text removed (it was faux-bolding the glyphs)
  * 'Orcestrator' typo
  * line positions of multi-line text recomputed from the new size. Many lines
    were stored at x=0 y=0 (Inkscape re-lays them out on load, browsers and PDF
    exporters do not), so they are rewritten explicitly
  * the per-element geometry tweaks in GEOMETRY / TEXT_MOVES, where larger text
    would otherwise overflow its box
  * v2 (2026-09-30): the Modes row is removed from the supervisor box (modes now
    get their own panel b, see make_fig1.py); the box is shortened to end level
    with the Browser box, its State & Coordination block re-centred, the dashed
    arrow to Tools lengthened to match, and the empty page margin under the
    drawing trimmed
  * v5 (2026-09-30): right-hand margins made equal to the left ones. The canvas is
    widened (VIEW_W) and the dashed container extended so Skills sits 3.95 units
    inside it and the container 4.22 units inside the page, exactly as Biologist
    and the container do on the left; the Watchdog is re-centred on the container.
    Font sizes keep their v4 size relative to the drawing (K_FONT), so every fit
    above still holds; printed sizes are 2.2 % smaller (smallest text 5.4 pt)
  * v10 (2026-10-01): the 14 icons traced from bitmaps of unknown origin are replaced
    by Tabler Icons (MIT, pinned in icons/tabler/, see ICON_SOURCES.md). Each new icon
    goes into the traced path's parent group, centred on and sized to the traced
    path's bounding box (measured in the source drawing). The Fiji / napari / Python
    logos (embedded PNGs) are untouched
"""
import io
import re
import xml.etree.ElementTree as ET

SRC, DST = "../Overview.svg", "Overview_nature.svg"
# Margins on the left: page -> container 4.22, container -> Biologist 3.95.
# Right side gets the same: Skills ends at 283.54 -> container 287.49 -> page 291.71.
LEFT_PAGE, LEFT_INNER = 4.22, 3.95
SKILLS_RIGHT = 242.41 + 41.13
CONTAINER_RIGHT = SKILLS_RIGHT + LEFT_INNER
VIEW_W = CONTAINER_RIGHT + LEFT_PAGE
K = 180.0 / VIEW_W        # mm per user unit at print size
K_FONT = 180.0 / 285.48001  # font calibration of v1-v4, keeps text/box proportions
VIEW_H = 184.0            # user units; the drawing ends at 181.4, the rest was white margin
LINE_HEIGHT = 1.25  # Inkscape default


def pt(p):
    return p * 0.3528 / K_FONT


# original size (user units) -> printed point size
SIZE_MAP = [
    (2.10, 5.5),   # subagent / librarian descriptions (printed 3.6 pt before)
    (2.75, 6.5),   # subagent titles, skill items, Librarian
    (2.90, 6.5),   # tool labels, mode pills, results, ledger, chat UI, recipe DB
    (3.10, 7.0),   # 'Provides input + task'
    (3.45, 8.0),   # Watchdog title
    (3.60, 7.0),   # State & coordination, Modes, Context line, legend
    (3.75, 9.0),   # Agentic-J title
    (4.00, 8.0),   # section headers
]

NS = {
    "": "http://www.w3.org/2000/svg",
    "inkscape": "http://www.inkscape.org/namespaces/inkscape",
    "sodipodi": "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd",
    "xlink": "http://www.w3.org/1999/xlink",
}
for k, v in NS.items():
    ET.register_namespace(k, v)
SVG = "{http://www.w3.org/2000/svg}"
ROLE = "{http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd}role"

# Geometry tweaks: rect id -> deltas (dx, dy, dw, dh) in user units, pre-shrink.
GEOMETRY = {
    "rect11": (0, 0, CONTAINER_RIGHT - (4.22 + 280.31), 0),  # dashed container, right edge
    "rect48-86-9-0": (-2.5, 0, 8.0, 0),      # Watchdog box: fit 'Controls Code Execution'
    "rect10-2": (-1.85, -0.5, 3.7, 1.6),     # Chat / User Interface box
    # supervisor box (inside g109, y scaled by 1.0976265): end level with the Browser
    # box bottom (106.05) instead of 118.68 now that the Modes row is gone
    "rect6": (0, 0, 0, (106.05 - 118.68324) / 1.0976265),
    "rect110": (0, 2.0 / 1.0976265, 0, 0),  # State & Coordination block, re-centred
    "rect6-2-3": (0, 2.0, 0, 0),             # Results
    "rect6-2-3-0": (0, 2.0, 0, 0),           # State Ledger + Task
}
# Elements removed outright: the Modes label and the three mode pills (with labels).
REMOVE = {"text30-5-2", "g116", "g115", "g114"}
# Path 'd' overrides. Dashed supervisor -> Tools arrow now starts at the shorter box.
PATH_D = {"path186": "m 130.14979,106.05 0.0105,22.50242"}
# Text moves: text id -> (dx, dy) in the element's own coordinates.
TEXT_MOVES = {
    "text30": (0, 0.9),                      # Agentic-J title off the top border
    "text19-8-74-8": (0, 1.0),               # State Ledger + Task: centred, +2 with its box
    "text30-5": (0, 2.0),                    # 'State & Coordination', moves with its block
    "text19-8-74": (0, 2.0),                 # 'Results', moves with its box
    # subagent titles / descriptions: shift left so the larger text fits the box
    "text58": (-2.0, 0.5), "text61-1": (-2.0, 1.0),        # Visual Judge
    "text56": (-2.0, 0.5), "text61-4": (-2.0, 0.5),        # Quality Assurance
    "text52": (-2.0, 0.5), "text61": (-2.0, 1.0),          # Discovery
    "text53": (-3.5, 0.5), "text61-5": (-14.18, 0.5),      # Fiji Coder: align desc with title
    "text55": (-2.0, 0.5), "text61-2": (-2.0, 1.0),        # Python Coder
    "text19-8-7-4-25": (0, 1.2),             # Recipe & Mistakes DB, now two lines
    **{t: (-1.5, 0) for t in ["text52-2", "text52-2-2", "text52-2-0", "text52-2-3",
                              "text52-2-00", "text52-2-6"]},  # skill items
}
# Per-element point size, overriding SIZE_MAP.
SIZE_OVERRIDE = {t: 6.0 for t in ["text52-2", "text52-2-2", "text52-2-0", "text52-2-3",
                                  "text52-2-00", "text52-2-6"]}
# Centre-anchored text switched to left-anchored (so it can align with its title).
ANCHOR_START = {"text61-5"}
# Rewrapped text: text id -> lines. Wording unchanged, only fewer lines so it fits.
RELINE = {
    "text61-1": ["Inspects inputs,", "segmentations & plots"],
    "text61": ["Overview of all installed", "packages and plugins"],
    "text61-2": ["Python packages,", "statistics and plotting"],
    "text19-8-7-4-25": ["Recipe &", "Mistakes DB"],
}


def target_size(orig):
    for limit, p in SIZE_MAP:
        if orig <= limit:
            return pt(p)
    return pt(8.0)


def sizes_in(el):
    return [float(x) for x in re.findall(r"font-size:([0-9.]+)px", el.get("style", ""))]


def set_style(el, key, val):
    st = el.get("style", "")
    if re.search(rf"(^|;){key}:", st):
        st = re.sub(rf"(^|;){key}:[^;]*", rf"\g<1>{key}:{val}", st)
    else:
        st = (st + ";" if st else "") + f"{key}:{val}"
    el.set("style", st)


tree = ET.parse(SRC)
root = tree.getroot()
root.set("width", "180mm")
root.set("height", f"{VIEW_H * K:.2f}mm")
root.set("viewBox", f"0 0 {VIEW_W:.4f} {VIEW_H}")

# Empty <text> leftovers from Inkscape editing: invisible, but they carry odd fonts
# and sizes that a journal preflight would flag.
parents = {c: p for p in root.iter() for c in p}
for el in list(root.iter()):
    if el.get("id") in REMOVE:
        parents[el].remove(el)
for text in list(root.iter(SVG + "text")):
    if not "".join(text.itertext()).strip():
        parents[text].remove(text)

for text in root.iter(SVG + "text"):
    tid = text.get("id")
    all_sizes = [s for e in text.iter() for s in sizes_in(e)]
    base = max(all_sizes)
    new = pt(SIZE_OVERRIDE[tid]) if tid in SIZE_OVERRIDE else target_size(base)
    ratio = new / base
    for e in text.iter():
        st = e.get("style")
        if not st:
            continue
        st = re.sub(r"font-size:([0-9.]+)px", lambda m: f"font-size:{float(m.group(1)) * ratio:.4f}px", st)
        st = st.replace("font-family:sans-serif", "font-family:Arial").replace("font-family:sans", "font-family:Arial")
        st = re.sub(r"-inkscape-font-specification:'?sans( Bold)?'?",
                    lambda m: "-inkscape-font-specification:'Arial%s'" % (m.group(1) or ""), st)
        st = re.sub(r"(^|;)stroke:#[0-9a-fA-F]{6}", r"\1stroke:none", st)
        e.set("style", st)

    if tid in ANCHOR_START:
        for e in text.iter():
            if e.get("style"):
                e.set("style", e.get("style").replace("text-anchor:middle", "text-anchor:start")
                      .replace("text-align:center", "text-align:start"))

    if tid in RELINE:
        lines = [t for t in text if t.tag == SVG + "tspan"]
        proto = lines[0]
        for t in lines:
            text.remove(t)
        for i, s in enumerate(RELINE[tid]):
            t = ET.SubElement(text, SVG + "tspan", dict(proto.attrib))
            t.set("id", f"{tid}-l{i}")
            t.text = s

    for e in text.iter():
        if e.text and "Orcestrator" in e.text:
            e.text = e.text.replace("Orcestrator", "Orchestrator")

    dx, dy = TEXT_MOVES.get(tid, (0.0, 0.0))
    x0 = float(text.get("x")) + dx
    y0 = float(text.get("y")) + dy
    text.set("x", f"{x0:.4f}")
    text.set("y", f"{y0:.4f}")
    lines = [t for t in text if t.tag == SVG + "tspan" and t.get(ROLE) == "line"]
    for i, t in enumerate(lines):
        t.set("x", f"{x0:.4f}")
        t.set("y", f"{y0 + i * new * LINE_HEIGHT:.4f}")

for el in root.iter(SVG + "rect"):
    d = GEOMETRY.get(el.get("id"))
    if d:
        for key, delta in zip(("x", "y", "width", "height"), d):
            el.set(key, f"{float(el.get(key)) + delta:.4f}")

for el in root.iter(SVG + "path"):
    if el.get("id") in PATH_D:
        el.set("d", PATH_D[el.get("id")])

# Re-centre the Watchdog group (box, icon, texts) on the widened container.
WATCHDOG_CX = 114.99 - 2.5 + (49.08 + 8.0) / 2
CONTAINER_CX = (LEFT_PAGE + CONTAINER_RIGHT) / 2
for el in root.iter(SVG + "g"):
    if el.get("id") == "g206":
        el.set("transform", f"translate({CONTAINER_CX - WATCHDOG_CX:.4f},0)")

# ── icon replacement (v10) ───────────────────────────────────────────────────
# traced path id -> (Tabler icon, bbox x, y, w, h of the traced path in root user
# units of ../Overview.svg, measured with chromium getBoundingClientRect)
ICONS = {
    "path148": ("user-shield", 118.04, 5.23, 9.01, 9.01),               # Watchdog
    "path1": ("eye", 186.53, 35.13, 11.64, 5.78),                       # Visual Judge
    "path137": ("zoom-check", 187.34, 51.92, 10.02, 10.02),             # Quality Assurance
    "path138": ("brain", 188.27, 68.96, 8.16, 9.77),                    # Discovery
    "path145": ("chart-dots", 245.47, 70.96, 6.46, 6.46),               # Plotting guidelines
    "path146": ("chart-histogram", 245.19, 81.26, 7.02, 7.02),          # Statistics
    "path147": ("building-bank", 245.00, 92.12, 7.41, 7.43),            # Publication standards
    "path139": ("database", 190.10, 137.25, 9.10, 9.08),                # Recipe & Mistakes DB
    "path136": ("book", 188.31, 160.15, 10.56, 9.53),                   # Librarian
    "path149": ("folder-open", 87.07, 148.22, 8.63, 7.91),              # File tools
    "path150": ("device-desktop-analytics", 102.33, 147.48, 9.41, 9.40),  # Inspect UI
    "path152": ("list-check", 124.62, 147.94, 8.47, 8.47),               # Check environment
    "path153": ("database-cog", 146.20, 147.94, 8.34, 8.62),            # Collect metadata
    "path151": ("code", 162.39, 147.48, 10.25, 9.40),                   # Script execution
}
ICON_DIR = "icons/tabler"
ICON_CONTENT = 18.0       # Tabler draws within ~3..21 of its 24-unit grid
ICON_STROKE = 1.25        # in Tabler units (their default is 2); matches the traced weight


def _affine(tr):
    m = [1, 0, 0, 1, 0, 0]
    for name, args in re.findall(r"(\w+)\(([^)]*)\)", tr or ""):
        a = [float(v) for v in re.split(r"[ ,]+", args.strip())]
        if name == "translate":
            n = [1, 0, 0, 1, a[0], a[1] if len(a) > 1 else 0]
        elif name == "scale":
            n = [a[0], 0, 0, a[1] if len(a) > 1 else a[0], 0, 0]
        elif name == "matrix":
            n = a
        else:
            raise ValueError(f"unhandled transform {name}")
        m = [m[0] * n[0] + m[2] * n[1], m[1] * n[0] + m[3] * n[1],
             m[0] * n[2] + m[2] * n[3], m[1] * n[2] + m[3] * n[3],
             m[0] * n[4] + m[2] * n[5] + m[4], m[1] * n[4] + m[3] * n[5] + m[5]]
    return m


def _tabler_children(name):
    """Drawing elements of a Tabler SVG (comment header and invisible frame dropped)."""
    src = io.open(f"{ICON_DIR}/{name}.svg", encoding="utf-8").read()
    src = re.sub(r"<!--.*?-->", "", src, flags=re.S)
    icon = ET.fromstring(src)
    return [c for c in icon if not (c.get("stroke") == "none" and c.get("d", "").startswith("M0 0h24"))]


parents = {c: p for p in root.iter() for c in p}
for el in list(root.iter(SVG + "path")):
    spec = ICONS.get(el.get("id"))
    if not spec:
        continue
    name, bx, by, bw, bh = spec
    par = parents[el]
    # parent's cumulative transform (root <- parent); the root svg itself has none
    chain, p = [], par
    while p is not None and p is not root:
        chain.append(p.get("transform"))
        p = parents.get(p)
    m = [1, 0, 0, 1, 0, 0]
    for tr in reversed(chain):
        n = _affine(tr)
        m = [m[0] * n[0] + m[2] * n[1], m[1] * n[0] + m[3] * n[1],
             m[0] * n[2] + m[2] * n[3], m[1] * n[2] + m[3] * n[3],
             m[0] * n[4] + m[2] * n[5] + m[4], m[1] * n[4] + m[3] * n[5] + m[5]]
    if abs(m[1]) > 1e-9 or abs(m[2]) > 1e-9:
        raise ValueError(f"{el.get('id')}: rotated parent not supported")
    cx, cy = bx + bw / 2, by + bh / 2                      # root coords
    lx, ly = (cx - m[4]) / m[0], (cy - m[5]) / m[3]        # parent-local coords
    sc = max(bw, bh) / ICON_CONTENT
    g = ET.Element(SVG + "g", {
        "id": f"icon-{name}",
        "transform": f"translate({lx - 12 * sc / m[0]:.4f},{ly - 12 * sc / m[3]:.4f}) "
                     f"scale({sc / m[0]:.5f},{sc / m[3]:.5f})",
        "fill": "none", "stroke": "#000000", "stroke-width": str(ICON_STROKE),
        "stroke-linecap": "round", "stroke-linejoin": "round",
    })
    for c in _tabler_children(name):
        for k in ("stroke", "fill"):
            if c.get(k) == "currentColor":
                del c.attrib[k]
        g.append(c)
    idx = list(par).index(el)
    par.remove(el)
    par.insert(idx, g)

tree.write(DST, encoding="UTF-8", xml_declaration=True)
print(f"wrote {DST}  (K={K:.4f})")
