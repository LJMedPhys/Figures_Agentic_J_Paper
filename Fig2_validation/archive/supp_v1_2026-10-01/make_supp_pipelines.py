"""Supplementary figures (v1, 2026-10-01): the step-by-step pipeline tables that Fig. 2 replaced with its scorecard.

  supplementary/Supp_CD177_pipelines.{pdf,svg,png}      = panel a of ../CD177_draft/make_figure.py (v7)
  supplementary/Supp_Parasite_pipelines.{pdf,svg,png}   = panel a of ../Parasite_draft/make_figure.py

Step texts, times and CD177 verdicts are copied verbatim from the drafts. Parasite verdict numbers are not
re-implemented here: the parasite data block of make_fig2.py is executed as is, so the supplement and Fig. 2
always show the same values. Verdict wording matches the Fig. 2 scorecard (slightly shortened from the draft
so the table fits 180 mm). Style as Fig. 2: Nature sizes (min 5 pt), 180 mm wide. Mounts: see run_supp.sh."""
import os, re
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


def mark(a, x, y, key, ok):
    a.text(x, y, '✓' if ok == 'ok' else '✗', ha='center', va='center', fontsize=7, fontweight='bold',
           color=COL[key] if ok == 'ok' else BAD, fontfamily='DejaVu Sans')


def legend_row(fig, a, x, y, items, extra=True):
    rend = fig.canvas.get_renderer()
    for k, t in items:
        a.add_patch(Rectangle((x, y - 1.1), 2.4, 2.2, fc=FILL[k], ec='k', lw=0.3)); tx = a.text(x + 3.1, y, t, va='center', fontsize=FS)
        x += 3.1 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0
    if extra:
        for ok, t in (('ok', 'published conclusion reproduced'), ('x', 'not reproduced')):
            a.text(x + 1, y, '✓' if ok == 'ok' else '✗', ha='center', va='center', fontsize=7, fontweight='bold',
                   color='#444444' if ok == 'ok' else BAD, fontfamily='DejaVu Sans')
            tx = a.text(x + 3, y, t, va='center', fontsize=FS); x += 3 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0


def cells(a, x0, cw, gap, y, rh, cl):
    for i, (tt, kind) in enumerate(cl):
        xx = x0 + i * (cw + gap)
        a.add_patch(Rectangle((xx, y), cw, rh, fc=FILL[kind], ec='k', lw=0.35))
        a.text(xx + cw / 2, y + rh / 2, tt, ha='center', va='center', fontsize=FS, linespacing=1.0)


# =============================== CD177 ===============================
H = 47; fig, ax_mm = new_fig(H)
a = ax_mm(0, 0, W, H); a.axis('off'); a.set_xlim(0, W); a.set_ylim(H, 0)
steps = ['Cell\nsegmentation', 'Nucleus\nsegmentation', 'QC &\npost-processing', 'Features', 'Statistics']
rows = [('Published', 'published', [('Cellpose (default)', 'manual'), ('ilastik pixel\nclassifier', 'manual'), ('dilation, crop,\nsize filter', 'manual'), ('25 hand-picked', 'manual'), ('cell-level\nMann–Whitney', 'manual')]),
        ('Codex CLI', 'codex', [('Cellpose cyto3\n(half resolution)', 'other'), ('erosion/threshold\nwithin cell', 'other'), ('shape QC', 'other'), ('30 features', 'other'), ('donor-level\nWilcoxon, FDR', 'other')]),
        ('Exploratory\nAgentic-J', 'fast', [('Cellpose cyto3', 'auto'), ('Otsu threshold\nwithin cell', 'auto'), ('border, size,\nsolidity filter', 'auto'), ('25 features', 'auto'), ('donor-level\nWilcoxon, FDR', 'auto')]),
        ('Conventional\nAgentic-J', 'automated', [('Cellpose cyto3\non optical density', 'auto'), ('threshold\nwithin cell', 'auto'), ('size/shape QC', 'auto'), ('32 features', 'auto'), ('donor-level\npermutation, FDR', 'auto')]),
        ('Conventional\nAgentic-J\n+ human', 'hil', [('Cellpose-SAM\n(user-specified)', 'hf'), ('micro-sam fine-tuned on\n11 hand-annotated tiles', 'hf'), ('lobe-to-cell\nassignment', 'auto'), ('11 features', 'auto'), ('image-level\nMann–Whitney, FDR', 'auto')])]
TIME = {'published': (3 * 7 * 24 * 60, 'weeks'), 'hil': (488, '488 min'), 'automated': (67, '67 min'), 'fast': (35, '35 min'), 'codex': (13, '13 min')}
VERD = {'published': ('ref', 'no meaningful difference'), 'hil': ('ok', '0/11 features significant'),
        'automated': ('ok', '0/32 significant (3/3 runs)'), 'fast': ('x', '5/25 significant (3/3 runs)'), 'codex': ('ok', '0/30 features significant')}
x0, cw, gap, rh, y0 = 13, 19.4, 0.5, 5.8, 6.0
STEP_END = x0 + 5 * (cw + gap) - gap
for i, s_ in enumerate(steps): a.text(x0 + i * (cw + gap) + cw / 2, 3.1, s_, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
TX, VX = STEP_END + 1.2, STEP_END + 33.5          # time bars, verdict column
a.text(TX + 14, 3.1, 'Time to result', ha='center', va='center', fontsize=FS, fontweight='bold')
a.text(VX + 16, 3.1, 'Scientific verdict', ha='center', va='center', fontsize=FS, fontweight='bold')
for r, (name, key, cl) in enumerate(rows):
    y = y0 + r * (rh + 0.5)
    a.text(0.3, y + rh / 2, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=1.0)
    cells(a, x0, cw, gap, y, rh, cl)
    ok, txt = VERD[key]
    if ok != 'ref': mark(a, VX + 1, y + rh / 2, key, ok)
    a.text(VX + 3, y + rh / 2, txt, va='center', fontsize=FS, color=COL['published'] if ok == 'ref' else 'k')
ta = ax_mm(TX, y0, 28, 5 * (rh + 0.5) - 0.5)
for r, (name, key, _) in enumerate(rows):
    v, lab = TIME[key]; ta.barh(r, v, height=0.55, color=COL[key], lw=0); ta.text(v * 1.3, r, lab, va='center', fontsize=FS)
ta.set_xscale('log'); ta.set_xlim(5, 4e5); ta.set_ylim(4.5, -0.5); ta.set_yticks([]); ta.minorticks_off()
ta.set_xticks([10, 60, 1440, 10080]); ta.set_xticklabels(['10 min', '1 h', '1 d', '1 wk']); ta.spines['left'].set_visible(False)
legend_row(fig, a, x0, 43.6, [('manual', 'Manual (analyst)'), ('hf', 'Human feedback (user)'), ('auto', 'Automatic (Agentic-J)'), ('other', 'Automatic (other agent)')])
for ext in ['pdf', 'svg', 'png']: fig.savefig(f'{OUT}/Supp_CD177_pipelines.{ext}', dpi=450 if ext == 'png' else None)
plt.close(fig)

# =============================== Parasite ===============================
# run the parasite data block of make_fig2.py verbatim (single source of truth for the numbers)
src = open('/work/make_fig2.py', encoding='utf-8').read()
block = src[src.index('P = \'/p\'; sys.path'):src.index('FS, FL = 5, 6')] + \
        src[src.index('# ======================= parasite data'):src.index('# ======================= a  scorecard')]
ns = {}
exec('import os, sys, csv, statistics as st\nimport numpy as np\nimport tifffile as tf\nfrom scipy.stats import mannwhitneyu\n' + block, ns)
agree, r1, S2, p2, K2, PUB, TOL = (ns[k] for k in ('agree', 'r1', 'S2', 'p2', 'K2', 'PUB', 'TOL'))

H = 48; fig, ax_mm = new_fig(H)
a = ax_mm(0, 0, W, H); a.axis('off'); a.set_xlim(0, W); a.set_ylim(H, 0)
steps = ['Registration\n(live → fixed)', 'Host / nucleus\nsegmentation', 'TUNEL call', 'Parasite\nlocalisation', 'Measurement\n& statistics']
rows = [('Published', 'published', [('by hand', 'manual'), ('hand-drawn\nROIs', 'manual'), ('by inspection', 'manual'), ('15 ROIs + 3\nbackgrounds\neach, by hand', 'manual'), ('ROI Manager →\nExcel;\nMann–Whitney', 'manual')]),
        ('Codex CLI', 'codex', [('SIFT + RANSAC,\nECC refined\n(1.6 px)', 'other'), ('T1 DAPI\nthreshold; T2\nCellpose cyto3', 'other'), ('largest-gap\ncutoff', 'other'), ('red-channel\npeaks + water-\nshed: 10 kept', 'other'), ('annulus\nbackground;\nno test', 'other')]),
        ('Exploratory\nAgentic-J', 'fast', [('TurboReg rigid\non brightfield', 'auto'), ('T1 StarDist;\nT2 Cellpose\ncyto3 (7 hosts)', 'auto'), ('largest-gap\ncutoff', 'auto'), ('C2 + C3 spot\nevidence:\n22 kept', 'auto'), ('annulus\nbackground;\nno test', 'auto')]),
        ('Conventional\nAgentic-J', 'automated', [('TurboReg /\nStackReg rigid', 'auto'), ('T1 DoG + water-\nshed; T2 Cellpose\ncyto3 (6 hosts)', 'auto'), ('2-component\nGMM', 'auto'), ('micro-sam\nvit_b_lm:\n105 → 5 kept', 'auto'), ('annulus\nbackground;\nno test', 'auto')]),
        ('Conventional\nAgentic-J\n+ human', 'hil', [('lab ROIs given,\norigin\ntranslated', 'hf'), ('lab ROIs\n(given)', 'hf'), ('largest-gap\ncutoff', 'hf'), ('lab ROIs paired\n1:1 with\nbackgrounds', 'hf'), ('native 16-bit;\nMann–Whitney', 'hf')])]
TIME = {'published': 'not recorded', 'hil': '32 min / n.r.', 'automated': '55 / 85 min', 'fast': '22 / 69 min', 'codex': '18 min (both)'}
def pfmt(p): return f'p = {p:.2g}'
w2n = sum(e <= TOL for e in S2['hil run 1']['err']); sc = S2['codex']; pc = p2[K2.index('Codex')]
VER = {'published': [('ref', '2/6 nuclei TUNEL+'), ('ref', f'8 TUNEL+ vs 7 TUNEL−, {pfmt(p2[K2.index(PUB)])}')],
       'hil': [('ok', f'{agree["hil"]}/6 calls, r = {r1["hil"]:.2f}'), ('ok', f'15/15; {w2n}/15 within {TOL}; {pfmt(p2[K2.index("Agent + human (run 1)")])}')],
       'automated': [('ok', f'{agree["automated"]}/6 calls, r = {r1["automated"]:.2f}'), ('x', f'{S2["automated"]["found"]}/15 found (+{S2["automated"]["extra"]}), no TUNEL+')],
       'fast': [('ok', f'{agree["fast"]}/6 calls, r = {r1["fast"]:.2f}, 1 split'), ('x', f'{S2["fast"]["found"]}/15 found (+{S2["fast"]["extra"]}), no host')],
       'codex': [('ok', f'{agree["codex"]}/6 calls, r = {r1["codex"]:.2f}'),
                 ('ok' if sc['grp'] == sc['found'] else 'x', f'{sc["found"]}/15 found, hosts right' + (f'; {pfmt(pc)}' if pc else ''))]}
x0, cw, gap, rh, y0 = 14.5, 15.0, 0.5, 6.4, 6.4
STEP_END = x0 + 5 * (cw + gap) - gap
for i, s_ in enumerate(steps): a.text(x0 + i * (cw + gap) + cw / 2, 3.2, s_, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
xt = STEP_END + 1.0; vx = [xt + 15.0, xt + 15.0 + 33.5]
a.text(xt + 7, 3.2, 'Time\nT1 / T2', ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
for j, t in enumerate(['Task 1  TUNEL calls\n(6 nuclei, fixed frame)', 'Task 2  published Fig. 4C\n(15 parasites, last live frame)']):
    a.text(vx[j] + 16, 3.2, t, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
for r, (name, key, cl) in enumerate(rows):
    y = y0 + r * (rh + 0.5)
    a.text(0.3, y + rh / 2, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=1.0)
    cells(a, x0, cw, gap, y, rh, cl)
    a.text(xt + 7, y + rh / 2, TIME[key], ha='center', va='center', fontsize=FS, color='#555555' if key == 'published' else 'k')
    for j, (ok, txt) in enumerate(VER[key]):
        if ok != 'ref': mark(a, vx[j] + 1, y + rh / 2, key, ok)
        a.text(vx[j] + 3, y + rh / 2, txt, va='center', fontsize=FS, color=COL['published'] if ok == 'ref' else 'k')
legend_row(fig, a, x0, 44.8, [('manual', 'Manual'), ('hf', 'Human feedback'), ('auto', 'Automatic (Agentic-J)'), ('other', 'Automatic (other agent)')])
for ext in ['pdf', 'svg', 'png']: fig.savefig(f'{OUT}/Supp_Parasite_pipelines.{ext}', dpi=450 if ext == 'png' else None)
print('parasite: calls agree', agree, '| task2', {m: (v['found'], v['extra']) for m, v in S2.items()})
