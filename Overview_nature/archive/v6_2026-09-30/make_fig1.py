"""Compose Fig. 1: (a) the overview from fix_overview.py, (b) the Agentic-J mode tree.

Run fix_overview.py first (render.sh does both). Units of the outer SVG are mm at
final print size (180 mm wide). Panel b is drawn in the overview's visual language:
same pastel fills (opacity 0.387 or the solid pill colours), black 0.265-unit
outlines scaled to print, Arial, filled-triangle arrowheads.

Mode tree (as specified by Lukas, 2026-09-30, revised same day):
    Agentic-J -> Analysis  -> Script-based -> Conventional -> optional human feedback
                                           -> Exploratory  -> optional human feedback
                           -> UI mode       -> human does the steps (agent gives click
                                              instructions, guides, explains parameters)
              -> Education                  -> human does the steps
Every branch ends at the human (right-hand column). The two automatic leaves merge
into one dashed arrow (optional feedback); UI mode and Education merge into one solid
arrow (the human performs the steps).
Subtitles come from the code (src/imagentj/prompts.py): Conventional = the former
'advanced' mode (phased pipeline: plan -> process -> stats -> plots -> docs -> QA);
Exploratory = the former 'quick' mode (one self-contained operation, no planning,
statistics or QA); Education = the tutor mode (concepts + short live demos).
"""
import io
import xml.etree.ElementTree as ET

import fix_overview as fo  # re-uses namespace registration and constants

W = 180.0
A_H = fo.VIEW_H * fo.K          # height of panel a in mm
B_TOP = A_H + 2.5               # panel b starts here
H = B_TOP + 46.0
PT = 0.3528
SW = 0.265 * fo.K               # outline width matching the overview's boxes
INK = "#000000"
TXT = "#030204"
GREY = "#4d4d4d"

# overview palette
GREEN = ("#5fc666", 0.387)      # supervisor green
TEAL = ("#a3ddd1", 1.0)         # old 'Advanced' pill
TEAL_L = ("#a3ddd1", 0.5)
YELLOW = ("#e2e790", 1.0)       # old 'Fast' pill
BLUE = ("#a3afdd", 1.0)         # old 'Education' pill
WHITE = ("#ffffff", 1.0)

out = []


def box(x, y, w, h, fill, rx=1.6, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" '
               f'fill="{fill[0]}" fill-opacity="{fill[1]}" stroke="{INK}" stroke-width="{SW:.3f}"'
               f' stroke-linejoin="round"{d}/>')


def text(x, y, s, pt=6.5, weight="normal", anchor="middle", color=TXT, style="normal"):
    s = s.replace("&", "&amp;")
    out.append(f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial" font-size="{pt * PT:.3f}" '
               f'font-weight="{weight}" font-style="{style}" text-anchor="{anchor}" fill="{color}">{s}</text>')


def path(pts, arrow=True, dash=None, both=False, sw=0.3):
    d = "M " + " L ".join(f"{a:.2f} {b:.2f}" for a, b in pts)
    m = ' marker-end="url(#fig1-ah)"' if arrow else ""
    m += ' marker-start="url(#fig1-ah)"' if both else ""
    da = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(f'<path d="{d}" fill="none" stroke="{INK}" stroke-width="{sw}"{m}{da} stroke-linejoin="round"/>')


def node(x, cy, w, h, fill, title, sub=None, tpt=7.0):
    box(x, cy - h / 2, w, h, fill)
    if isinstance(sub, (list, tuple)):          # multi-line subtitle, block centred
        top = cy - (len(sub) * 2.55) / 2 + 0.9
        text(x + w / 2, top, title, tpt, "bold")
        for i, line in enumerate(sub):
            text(x + w / 2, top + 2.95 + i * 2.55, line, 6.0, color=GREY)
    elif sub:
        text(x + w / 2, cy - 0.35, title, tpt, "bold")
        text(x + w / 2, cy + 2.6, sub, 6.0, color=GREY)
    else:
        text(x + w / 2, cy + tpt * PT * 0.36, title, tpt, "bold")


def elbow(x0, y0, x1, y1):
    """Tree connector: parent right edge -> child left edge, orthogonal."""
    xm = (x0 + x1) / 2
    path([(x0, y0), (xm, y0), (xm, y1), (x1 - 0.2, y1)])


def human_end(starts, px, doer, xm):
    """Merge several leaves into ONE arrow to ONE person icon (mirror of the tree fork).

    starts: [(x_right_edge, cy), ...]; the branches join at x = xm, and a single
    arrow runs from the join to the person, centred between the rows."""
    dash = None if doer else "0.9 0.6"
    ys = [cy for _, cy in starts]
    ym = (min(ys) + max(ys)) / 2
    for x0, cy in starts:
        path([(x0, cy), (xm, cy)], arrow=False, dash=dash)
    path([(xm, min(ys)), (xm, max(ys))], arrow=False, dash=dash)
    path([(xm, ym), (px - 3.6, ym)], dash=dash)
    l1, l2 = ("human does", "the steps") if doer else ("optional", "feedback")
    person(px, ym - 0.4)
    text(px + 4.0, ym - 0.35, l1, 6.0, anchor="start", color=GREY)
    text(px + 4.0, ym + 2.25, l2, 6.0, anchor="start", color=GREY)


def person(cx, cy, s=1.0):
    out.append(f'<g fill="none" stroke="{INK}" stroke-width="0.3">'
               f'<circle cx="{cx:.2f}" cy="{cy - 1.9 * s:.2f}" r="{1.5 * s:.2f}"/>'
               f'<path d="M {cx - 2.9 * s:.2f} {cy + 3.1 * s:.2f} Q {cx - 2.9 * s:.2f} {cy + 0.2 * s:.2f} '
               f'{cx:.2f} {cy + 0.2 * s:.2f} Q {cx + 2.9 * s:.2f} {cy + 0.2 * s:.2f} '
               f'{cx + 2.9 * s:.2f} {cy + 3.1 * s:.2f}"/></g>')


# ── panel b: mode tree ────────────────────────────────────────────────────────
NH = 8.0                                   # node height
# row centres: the two automatic leaves carry 2-line subtitles, so they get 12 mm
r = [B_TOP + 6.8, B_TOP + 18.8, B_TOP + 30.0, B_TOP + 40.5]
y_script = (r[0] + r[1]) / 2
y_analysis = (y_script + r[2]) / 2
y_root = (y_analysis + r[3]) / 2

C0, W0 = 3.0, 24.0     # root
C1, W1 = 33.0, 31.0    # analysis / education
C2, W2 = 71.0, 36.0    # script-based / UI mode
C3, W3 = 114.0, 36.0   # conventional / exploratory
PX = 161.5             # human column; labels end inside panel a's container edge (177.4 mm)

node(C0, y_root, W0, 11.0, GREEN, "Agentic-J", "mode set per request", tpt=8.0)
node(C1, y_analysis, W1, NH, TEAL, "Analysis", "process & measure images")
node(C1, r[3], W1, NH, BLUE, "Education", "tutor: concepts & demos")
node(C2, y_script, W2, NH, TEAL_L, "Script-based", "agent writes & runs scripts")
node(C2, r[2], W2, 10.6, YELLOW, "UI mode",
     ["gives click instructions,", "guides & explains parameters"])
node(C3, r[0], W3, 10.6, WHITE, "Conventional",
     ["full pipeline: plan, stats, QA", "knowledge base, concept library"])
node(C3, r[1], W3, 10.6, WHITE, "Exploratory", ["minimal steps,", "lightweight harness"])

elbow(C0 + W0, y_root, C1, y_analysis)
elbow(C0 + W0, y_root, C1, r[3])
elbow(C1 + W1, y_analysis, C2, y_script)
elbow(C1 + W1, y_analysis, C2, r[2])
elbow(C2 + W2, y_script, C3, r[0])
elbow(C2 + W2, y_script, C3, r[1])

human_end([(C3 + W3, r[0]), (C3 + W3, r[1])], PX, doer=False, xm=C3 + W3 + 2.5)
human_end([(C2 + W2, r[2]), (C1 + W1, r[3])], PX, doer=True, xm=C3 + W3 + 2.5)

# ── assemble ──────────────────────────────────────────────────────────────────
ET.register_namespace("", "http://www.w3.org/2000/svg")
a_root = ET.parse("Overview_nature.svg").getroot()
for k in ("width", "height"):
    a_root.attrib.pop(k, None)
a_root.set("x", "0")
a_root.set("y", "0")
a_root.set("width", f"{W}")
a_root.set("height", f"{A_H:.3f}")
panel_a = ET.tostring(a_root, encoding="unicode")

svg = (f'<?xml version="1.0" encoding="UTF-8"?>\n'
       f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
       f'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
       f'xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd" '
       f'width="{W}mm" height="{H:.2f}mm" viewBox="0 0 {W} {H:.2f}">\n'
       '<defs><marker id="fig1-ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="4.4" '
       'markerHeight="4.4" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="#000"/></marker></defs>\n'
       f'<rect width="{W}" height="{H:.2f}" fill="#ffffff"/>\n'
       + panel_a + "\n"
       + f'<text x="1.2" y="4.6" font-family="Arial" font-size="{8 * PT:.3f}" font-weight="bold" fill="{TXT}">a</text>\n'
       + f'<text x="1.2" y="{B_TOP + 2.4:.2f}" font-family="Arial" font-size="{8 * PT:.3f}" font-weight="bold" fill="{TXT}">b</text>\n'
       + "\n".join(out) + "\n</svg>\n")
io.open("Fig1_overview_modes.svg", "w", encoding="utf-8").write(svg)
print(f"wrote Fig1_overview_modes.svg  ({W} x {H:.1f} mm)")
