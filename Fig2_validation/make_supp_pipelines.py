"""Supplementary figures (v1, 2026-10-01): the step-by-step pipeline tables that Fig. 2 replaced with its scorecard.

  supplementary/Supp_CD177_pipelines.{pdf,svg,png}      = panel a of ../CD177_draft/make_figure.py (v7)
  supplementary/Supp_Parasite_pipelines.{pdf,svg,png}   = panel a of ../Parasite_draft/make_figure.py

Step texts are copied verbatim from the drafts. v2: everything already shown in Fig. 2 (time to result,
verdicts, task 1 / task 2 outcomes) is cut, so these are pure method tables and need no run data; the five
step columns use the full width. Style as Fig. 2: Nature sizes (min 5 pt), 180 mm wide."""
import os, re, textwrap
import matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

FS, FL = 5, 6
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': FS,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'xtick.labelsize': FS, 'ytick.labelsize': FS, 'axes.labelsize': FS, 'axes.titlesize': FS, 'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5})
MM = 1 / 25.4; W = 180
COL = {'published': '#7b3294', 'hil': '#0571b0', 'automated': '#1a9641', 'fast': '#a6611a', 'codex': '#636363'}
FILL = {'manual': '#f7c873', 'hf': '#8ecdf2', 'auto': '#a6dba0', 'other': '#d9d9d9'}
BAD = '#c0392b'
OUT = '/work/supplementary'; os.makedirs(OUT, exist_ok=True)


def new_fig(h):
    fig = plt.figure(figsize=(W * MM, h * MM))
    ax_mm = lambda x, y, w, hh, **k: fig.add_axes([x / W, 1 - (y + hh) / h, w / W, hh / h], **k)
    return fig, ax_mm


def legend_row(fig, a, x, y, items, extra=False):
    rend = fig.canvas.get_renderer()
    for k, t in items:
        a.add_patch(Rectangle((x, y - 1.1), 2.4, 2.2, fc=FILL[k], ec='k', lw=0.3)); tx = a.text(x + 3.1, y, t, va='center', fontsize=FS)
        x += 3.1 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0
    if extra:
        for ok, t in (('ok', 'published conclusion reproduced'), ('x', 'not reproduced')):
            a.text(x + 1, y, '✓' if ok == 'ok' else '✗', ha='center', va='center', fontsize=7, fontweight='bold',
                   color='#444444' if ok == 'ok' else BAD, fontfamily='DejaVu Sans')
            tx = a.text(x + 3, y, t, va='center', fontsize=FS); x += 3 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0


def rewrap(tt, width=34):
    """The drafts broke lines for ~15 mm cells; re-break for the ~33 mm cells here (words unchanged):
    one line if it fits; else keep the draft's own break if it was already two lines; else break after
    each ';' (separates the T1 / T2 parts) and wrap. A break after a hyphen joined a split word
    ('water-\nshed' -> 'watershed')."""
    flat = re.sub(r'(?<=[a-z])-\n(?=[a-z])', '', tt).replace('\n', ' ')
    if len(flat) <= width: return flat
    if tt.count('\n') == 1 and '-\n' not in tt and all(len(l) <= width for l in tt.split('\n')): return tt
    return '\n'.join(balanced(part, width) for part in re.split(r'(?<=;) ', flat))


def balanced(text, width):
    """Two lines split as evenly as possible (no orphan word like '(1.6 / px)'); textwrap beyond two lines."""
    if len(text) <= width: return text
    w = text.split(' ')
    best = min((max(len(' '.join(w[:i])), len(' '.join(w[i:]))), i) for i in range(1, len(w)))
    return ' '.join(w[:best[1]]) + '\n' + ' '.join(w[best[1]:]) if best[0] <= width else textwrap.fill(text, width)


def cells(a, x0, cw, gap, y, rh, cl):
    for i, (tt, kind) in enumerate(cl):
        tt = rewrap(tt)
        xx = x0 + i * (cw + gap)
        a.add_patch(Rectangle((xx, y), cw, rh, fc=FILL[kind], ec='k', lw=0.35))
        a.text(xx + cw / 2, y + rh / 2, tt, ha='center', va='center', fontsize=FS, linespacing=1.0)


# =============================== CD177 ===============================
H = 43; fig, ax_mm = new_fig(H)
a = ax_mm(0, 0, W, H); a.axis('off'); a.set_xlim(0, W); a.set_ylim(H, 0)
steps = ['Cell\nsegmentation', 'Nucleus\nsegmentation', 'QC &\npost-processing', 'Features', 'Statistics']
rows = [('Published', 'published', [('Cellpose (default)', 'manual'), ('ilastik pixel\nclassifier', 'manual'), ('dilation, crop,\nsize filter', 'manual'), ('25 hand-picked', 'manual'), ('cell-level\nMann–Whitney', 'manual')]),
        ('Codex CLI', 'codex', [('Cellpose cyto3\n(half resolution)', 'other'), ('erosion/threshold\nwithin cell', 'other'), ('shape QC', 'other'), ('30 features', 'other'), ('donor-level\nWilcoxon, FDR', 'other')]),
        ('Exploratory\nAgentic-J', 'fast', [('Cellpose cyto3', 'auto'), ('Otsu threshold\nwithin cell', 'auto'), ('border, size,\nsolidity filter', 'auto'), ('25 features', 'auto'), ('donor-level\nWilcoxon, FDR', 'auto')]),
        ('Conventional\nAgentic-J', 'automated', [('Cellpose cyto3\non optical density', 'auto'), ('threshold\nwithin cell', 'auto'), ('size/shape QC', 'auto'), ('32 features', 'auto'), ('donor-level\npermutation, FDR', 'auto')]),
        ('Conventional\nAgentic-J\n+ human', 'hil', [('Cellpose-SAM\n(user-specified)', 'hf'), ('micro-sam fine-tuned on\n11 hand-annotated tiles', 'hf'), ('lobe-to-cell\nassignment', 'auto'), ('11 features', 'auto'), ('image-level\nMann–Whitney, FDR', 'auto')])]
x0, gap, rh, y0 = 13, 0.6, 5.8, 6.0
cw = (W - 0.6 - x0 - 4 * gap) / 5                # five step columns over the full width (0.6 mm right margin)
for i, s_ in enumerate(steps): a.text(x0 + i * (cw + gap) + cw / 2, 3.1, s_, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
for r, (name, key, cl) in enumerate(rows):
    y = y0 + r * (rh + 0.5)
    a.text(0.3, y + rh / 2, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=1.0)
    cells(a, x0, cw, gap, y, rh, cl)
legend_row(fig, a, x0, 40.0, [('manual', 'Manual (analyst)'), ('hf', 'Human feedback (user)'), ('auto', 'Automatic (Agentic-J)'), ('other', 'Automatic (other agent)')])
for ext in ['pdf', 'svg', 'png']: fig.savefig(f'{OUT}/Supp_CD177_pipelines.{ext}', dpi=450 if ext == 'png' else None)
plt.close(fig)

# =============================== Parasite ===============================
H = 44; fig, ax_mm = new_fig(H)
a = ax_mm(0, 0, W, H); a.axis('off'); a.set_xlim(0, W); a.set_ylim(H, 0)
steps = ['Registration\n(live → fixed)', 'Host / nucleus\nsegmentation', 'TUNEL call', 'Parasite\nlocalisation', 'Measurement\n& statistics']
rows = [('Published', 'published', [('by hand', 'manual'), ('hand-drawn\nROIs', 'manual'), ('by inspection', 'manual'), ('15 ROIs + 3\nbackgrounds\neach, by hand', 'manual'), ('ROI Manager →\nExcel;\nMann–Whitney', 'manual')]),
        ('Codex CLI', 'codex', [('SIFT + RANSAC,\nECC refined\n(1.6 px)', 'other'), ('T1 DAPI\nthreshold; T2\nCellpose cyto3', 'other'), ('largest-gap\ncutoff', 'other'), ('red-channel\npeaks + water-\nshed: 10 kept', 'other'), ('annulus\nbackground;\nno test', 'other')]),
        ('Exploratory\nAgentic-J', 'fast', [('TurboReg rigid\non brightfield', 'auto'), ('T1 StarDist;\nT2 Cellpose\ncyto3 (7 hosts)', 'auto'), ('largest-gap\ncutoff', 'auto'), ('C2 + C3 spot\nevidence:\n22 kept', 'auto'), ('annulus\nbackground;\nno test', 'auto')]),
        ('Conventional\nAgentic-J', 'automated', [('TurboReg /\nStackReg rigid', 'auto'), ('T1 DoG + water-\nshed; T2 Cellpose\ncyto3 (6 hosts)', 'auto'), ('2-component\nGMM', 'auto'), ('micro-sam\nvit_b_lm:\n105 → 5 kept', 'auto'), ('annulus\nbackground;\nno test', 'auto')]),
        ('Conventional\nAgentic-J\n+ human', 'hil', [('lab ROIs given,\norigin\ntranslated', 'hf'), ('lab ROIs\n(given)', 'hf'), ('largest-gap\ncutoff', 'hf'), ('lab ROIs paired\n1:1 with\nbackgrounds', 'hf'), ('native 16-bit;\nMann–Whitney', 'hf')])]
x0, gap, rh, y0 = 14.5, 0.6, 5.8, 6.4
cw = (W - 0.6 - x0 - 4 * gap) / 5
for i, s_ in enumerate(steps): a.text(x0 + i * (cw + gap) + cw / 2, 3.2, s_, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
for r, (name, key, cl) in enumerate(rows):
    y = y0 + r * (rh + 0.5)
    a.text(0.3, y + rh / 2, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=1.0)
    cells(a, x0, cw, gap, y, rh, cl)
legend_row(fig, a, x0, 40.6, [('manual', 'Manual'), ('hf', 'Human feedback'), ('auto', 'Automatic (Agentic-J)'), ('other', 'Automatic (other agent)')])
for ext in ['pdf', 'svg', 'png']: fig.savefig(f'{OUT}/Supp_Parasite_pipelines.{ext}', dpi=450 if ext == 'png' else None)
print('wrote supplementary/Supp_CD177_pipelines + Supp_Parasite_pipelines')
