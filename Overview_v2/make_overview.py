"""Generate Overview_v2.svg: the Agentic-J architecture figure at Nature Methods size.

Edit this script, not the SVG. Units are millimetres (viewBox = mm), sized for a
180 mm double-column figure. Text stays >= 5.5 pt; body text is 6-6.5 pt Arial.
Render with ./render.sh (chromium from the presenton image, host Liberation Sans).
"""
import io

W, H = 180.0, 123.0
PT = 0.3528  # mm per point
FONT = "Arial, 'Liberation Sans', Helvetica, sans-serif"

# Muted palette, one hue per layer (no red/green pairing).
INK = "#1a1a1a"
GREY = "#6b7580"
USER = ("#eef1f4", "#55606b")
UI = ("#dde8f3", "#3e6a94")
SUP = ("#e8e2f3", "#5d4a91")
SUP_IN = "#f7f5fb"
LOOP = ("#fdf1dc", "#b8741a")
VERIFY = "#f7d9a6"
FAIL = "#c8531b"
KNOW = ("#e0eeeb", "#3b7f74")
LLM = ("#f1f1f1", "#7a7a7a")

out = []


def rect(x, y, w, h, fill, stroke, rx=1.6, sw=0.25, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" '
               f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')


def text(x, y, s, pt=6.0, weight="normal", anchor="start", color=INK, style="normal", rot=None):
    t = f' transform="rotate({rot} {x:.2f} {y:.2f})"' if rot is not None else ""
    s = s.replace("&", "&amp;")
    out.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{pt * PT:.3f}" font-weight="{weight}" '
               f'font-style="{style}" text-anchor="{anchor}" fill="{color}"{t}>{s}</text>')


def lines(x, y, rows, pt=6.0, lh=None, **kw):
    lh = lh or pt * PT * 1.25
    for i, r in enumerate(rows):
        text(x, y + i * lh, r, pt, **kw)


def arrow(pts, color=INK, sw=0.3, both=False, dash=None):
    d = "M " + " L ".join(f"{a:.2f} {b:.2f}" for a, b in pts)
    m = "ah-fail" if color == FAIL else "ah"
    ms = f' marker-start="url(#{m})"' if both else ""
    da = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw}" '
               f'marker-end="url(#{m})"{ms}{da} stroke-linejoin="round"/>')


def person(cx, cy, s=1.0, color=USER[1]):
    out.append(f'<g fill="none" stroke="{color}" stroke-width="0.35">'
               f'<circle cx="{cx}" cy="{cy - 2.2 * s}" r="{1.8 * s}"/>'
               f'<path d="M {cx - 3.4 * s} {cy + 3.6 * s} Q {cx - 3.4 * s} {cy + 0.2 * s} {cx} {cy + 0.2 * s} '
               f'Q {cx + 3.4 * s} {cy + 0.2 * s} {cx + 3.4 * s} {cy + 3.6 * s}"/></g>')


def step(n, x, y, w, h, title, rows, fill="#ffffff"):
    rect(x, y, w, h, fill, LOOP[1], rx=1.2, sw=0.22)
    out.append(f'<circle cx="{x + 3.6:.2f}" cy="{y + h / 2:.2f}" r="2.3" fill="{LOOP[1]}"/>')
    text(x + 3.6, y + h / 2 + 0.85, str(n), 6.5, "bold", "middle", "#ffffff")
    text(x + 7.2, y + 3.9, title, 6.5, "bold")
    lines(x + 7.2, y + 7.2, rows, 6.0, lh=2.75)


# ── outside the container: user ───────────────────────────────────────────────
rect(2, 38, 28, 36, *USER)
text(16, 43.5, "User", 8, "bold", "middle")
person(16, 51)
lines(16, 60.5, ["task and images,", "in plain language"], 6.0, anchor="middle")
lines(16, 67.5, ["no coding needed"], 6.0, anchor="middle", color=GREY, style="italic")

text(3, 81, "Human in the loop", 6.5, "bold")
lines(3, 85, ["• approves the pipeline", "• approves plugin installs", "• rates the results",
              "• annotates tiles for"], 6.0, lh=2.8)
text(5.2, 85 + 4 * 2.8, "micro-sam fine-tuning", 6.0)

# ── container boundary ────────────────────────────────────────────────────────
rect(33, 4, 144, 106, "none", GREY, rx=2.5, sw=0.3, dash="1.4 0.9")
text(35.5, 8.6, "Docker container  ·  sandboxed  ·  reproducible  ·  optional GPU", 6.0,
     color=GREY, style="italic")

# ── virtual desktop ───────────────────────────────────────────────────────────
rect(36, 26, 29, 56, *UI)
text(50.5, 31, "Virtual desktop", 7, "bold", "middle")
text(50.5, 34.6, "in the browser (noVNC)", 6.0, anchor="middle", color=GREY)
for i, (name, sub) in enumerate([("Chat panel", "talk to the agent"),
                                 ("Fiji", "agent's scripts run live"),
                                 ("napari", "viewer + plugins")]):
    y = 38 + i * 14.3
    rect(38.5, y, 24, 11.5, "#ffffff", UI[1], rx=1.1, sw=0.2)
    text(50.5, y + 4.9, name, 6.5, "bold", "middle")
    text(50.5, y + 8.5, sub, 5.5, anchor="middle", color=GREY)

# ── supervisor ────────────────────────────────────────────────────────────────
rect(70, 12, 46, 70, *SUP)
text(93, 18, "Supervisor", 8, "bold", "middle")
text(93, 22.2, "plans · delegates · integrates", 6.0, anchor="middle")
text(73, 27.2, "Mode", 6.0, color=GREY)
for i, m in enumerate(["Quick", "Advanced", "Education"]):
    x = 73 + i * 13.8
    rect(x, 28.5, 12.3, 5.6, "#ffffff", SUP[1], rx=2.8, sw=0.2)
    text(x + 6.15, 32.3, m, 6.0, anchor="middle")

rect(73, 37.5, 40, 17, SUP_IN, SUP[1], rx=1.1, sw=0.2)
text(75, 42, "State ledger", 6.5, "bold")
lines(75, 45.6, ["persistent project state:", "metadata · channels ·", "append-only step log"], 6.0, lh=2.75)

rect(73, 58, 40, 21, SUP_IN, SUP[1], rx=1.1, sw=0.2)
text(75, 62.5, "Tools", 6.5, "bold")
lines(75, 66.1, ["file & metadata inspection", "UI inspection · screenshots",
                 "script execution · env check", "napari MCP · web search"], 6.0, lh=2.75)

# ── closed-loop analysis ──────────────────────────────────────────────────────
rect(124, 12, 51, 70, *LOOP)
text(149.5, 18, "Closed-loop analysis", 8, "bold", "middle")
text(149.5, 22.2, "subagents; nothing counts until verified", 6.0, anchor="middle")
SX, SW, SH = 127, 38, 11.5
ys = [25, 39, 53, 67]
step(1, SX, ys[0], SW, SH, "Discover", ["Plugin manager picks a", "Fiji, Python or napari tool"])
step(2, SX, ys[1], SW, SH, "Write code", ["Fiji coder: macro / Groovy", "Python analyst: stats, plots"])
step(3, SX, ys[2], SW, SH, "Execute", ["in live Fiji or a conda env;", "Watchdog stops stuck runs"])
step(4, SX, ys[3], SW, SH, "Verify", ["Vision judge: masks, plots", "QA gate: community standards"],
     fill=VERIFY)
cx = SX + SW / 2
for a in range(3):
    arrow([(cx, ys[a] + SH), (cx, ys[a + 1] - 0.4)])
# fail -> debugger -> back to code
fy, ty = ys[3] + SH / 2, ys[1] + SH / 2
arrow([(SX + SW, fy), (170, fy), (170, ty), (SX + SW + 0.4, ty)], color=FAIL, sw=0.35)
text(172.6, (fy + ty) / 2, "fail → debugger fixes", 6.0, "bold", "middle", FAIL, rot=-90)

# ── supervisor <-> desktop / loop ─────────────────────────────────────────────
arrow([(30, 56), (35.6, 56)], both=True)                    # user <-> desktop
arrow([(65, 43.8), (69.6, 43.8)], both=True)                # chat <-> supervisor
arrow([(116, 30.8), (126.6, 30.8)])                         # delegate
text(120, 29.4, "tasks", 5.5, anchor="middle", color=GREY)
arrow([(SX, 72.8), (116.4, 72.8)])                          # verified results back
text(120.5, 71.4, "results", 5.5, anchor="middle", color=GREY)
text(120.5, 76.4, "(typed)", 5.5, anchor="middle", color=GREY)

# ── knowledge layer ───────────────────────────────────────────────────────────
rect(36, 89, 54, 18, *KNOW)
text(38.5, 94, "Knowledge base", 7, "bold")
lines(38.5, 98.2, ["documentation (retrieval) · concept library", "bioimage-analysis teaching material",
                   "used for planning and Education mode"], 6.0, lh=2.8)

rect(94, 89, 38, 18, *KNOW)
text(96.5, 94, "Skills", 7, "bold")
lines(96.5, 98.2, ["Fiji · Python · napari plugins", "plotting · statistics",
                   "publication standards"], 6.0, lh=2.8)

rect(136, 89, 39, 18, *KNOW)
text(138.5, 94, "Recipe & mistake memory", 7, "bold")
lines(138.5, 98.2, ["curated by the Librarian", "after verified runs;",
                    "recalled before coding"], 6.0, lh=2.8)

arrow([(63, 89), (63, 85), (80, 85), (80, 82.4)])           # knowledge -> supervisor
arrow([(106, 89), (106, 82.4)])                             # skills -> supervisor
arrow([(128.5, 89), (128.5, 82.4)])                         # skills -> loop
arrow([(144, 82), (144, 88.6)])                             # verified run -> memory
text(145.2, 86.3, "pass", 5.5, "bold")
arrow([(163, 89), (163, 82.4)])                             # memory -> loop
text(164.2, 86.3, "recall", 5.5, color=GREY)

# ── LLM backend (outside) ─────────────────────────────────────────────────────
rect(33, 113, 144, 7.5, *LLM)
text(36, 117.9, "LLM backend", 6.5, "bold")
text(56, 117.9, "one model per role, set in config (OpenRouter or OpenAI): supervisor · coders · "
               "analyst · vision judge · helpers", 6.0)
arrow([(105, 113), (105, 110.4)], color=GREY, dash="0.8 0.6")

svg = (f'<?xml version="1.0" encoding="UTF-8"?>\n'
       f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H}mm" '
       f'viewBox="0 0 {W} {H}" font-family="{FONT}">\n'
       '<defs>'
       '<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="4.2" markerHeight="4.2" '
       f'orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="{INK}"/></marker>'
       '<marker id="ah-fail" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="4" markerHeight="4" '
       f'orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="{FAIL}"/></marker>'
       '</defs>\n'
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>\n' + "\n".join(out) + "\n</svg>\n")
io.open("Overview_v2.svg", "w", encoding="utf-8").write(svg)
print("wrote Overview_v2.svg")
