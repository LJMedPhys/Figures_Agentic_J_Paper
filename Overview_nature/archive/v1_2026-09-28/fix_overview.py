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
"""
import io
import re
import xml.etree.ElementTree as ET

SRC, DST = "../Overview.svg", "Overview_nature.svg"
K = 180.0 / 285.48001
LINE_HEIGHT = 1.25  # Inkscape default


def pt(p):
    return p * 0.3528 / K


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
    "rect48-86-9-0": (-2.5, 0, 8.0, 0),      # Watchdog box: fit 'Controls Code Execution'
    "rect10-2": (-1.85, -0.5, 3.7, 1.6),     # Chat / User Interface box
    "rect6-2-3-1": (-0.75, 0, 1.5, 0),       # mode pills
    "rect6-2-3-1-6": (-0.75, 0, 1.5, 0),
    "rect6-2-3-1-2": (-0.75, 0, 1.5, 0),
}
# Text moves: text id -> (dx, dy) in the element's own coordinates.
TEXT_MOVES = {
    "text30": (0, 0.9),                      # Agentic-J title off the top border
    "text19-8-74-8": (0, -1.0),              # State Ledger + Task: centre in its box
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
root.set("height", f"{192.66 * K:.2f}mm")

# Empty <text> leftovers from Inkscape editing: invisible, but they carry odd fonts
# and sizes that a journal preflight would flag.
parents = {c: p for p in root.iter() for c in p}
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

tree.write(DST, encoding="UTF-8", xml_declaration=True)
print(f"wrote {DST}  (K={K:.4f})")
