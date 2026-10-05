"""Fig. 2, full-page variant (v1, 2026-10-02): BOTH original draft figures (CD177 + parasite) with every plot,
re-laid out on one page at journal size: 180 mm wide, <= 247 mm deep (legend on a separate page), all text >= 5 pt,
PNG at 600 dpi (combination artwork; PDF/SVG are vector).

  a  CD177 pipeline table (steps, time to result, verdict)          = CD177 draft a
  e  parasite pipeline table (steps, times, task 1 / task 2 verdict) = parasite draft a
  left column, CD177:  b gallery (4 cells), c agreement with published masks, d per-cell distributions
  right column, parasite: f gallery (field, cell 3, cell 6), g corrected TUNEL signal per nucleus,
                          h count comparison with the published analysis (= parasite draft c),
                          i proliferation index by host TUNEL status + delta median

Generated from make_fig2.py (same data code, same panels) by build_full.py-style replacements: the scorecard
is swapped for the two draft tables, the right column gets the count table, panels g and i are a little shorter.
Mounts: same as run.sh.
v2 (2026-10-02): re-arranged into two blocks, CD177 on top (a table; b gallery | c, d) and Leishmania below
(e table; f gallery + g | h + i), letters in reading order; Exploratory colour moved away from red; CD177 gallery
tiles re-rendered with a black 'published cell missed' outline (was red, too close to the orange 'not measured');
arrow on the Codex false nucleus fragment; human-in-the-loop parasite times 34 / 26 min; TUNEL+ / TUNEL- explained;
bracketed asides written out (axis units kept); task 1 verdict states which nuclei differ."""
import os, sys, csv, statistics as st
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, Ellipse, FancyBboxPatch
from matplotlib.lines import Line2D
import tifffile as tf
from skimage.measure import find_contours
from scipy.stats import mannwhitneyu

P = '/p'; sys.path.insert(0, f'{P}/scripts'); os.chdir(f'{P}/scripts')
import ij_roi, parasite_gt as gt
import figure_tasks23_dotplots as dp, score_auto_run as sa, build_workbook as bw
from figure_task1_overlays import load_run, stretch

FS, FL = 5, 6
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': FS,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'xtick.labelsize': FS, 'ytick.labelsize': FS, 'axes.labelsize': FS, 'axes.titlesize': FS, 'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5})
MM = 1 / 25.4; W, H = 180, 247
fig = plt.figure(figsize=(W * MM, H * MM))
def ax_mm(x, y, w, h, **k): return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **k)
def letter(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight='bold', va='top')
def head(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=FL, va='top')
COL = {'published': '#7b3294', 'hil': '#0571b0', 'automated': '#1a9641', 'fast': '#b8960c', 'codex': '#636363'}
FILL = {'manual': '#f7c873', 'hf': '#8ecdf2', 'auto': '#a6dba0', 'other': '#d9d9d9'}
SH = {'published': 'Published', 'hil': 'Conventional Agentic-J + human', 'automated': 'Conventional Agentic-J',
      'fast': 'Exploratory Agentic-J', 'codex': 'Codex'}
BAD = '#c0392b'
CD = '/cd177'

# ======================= parasite data (verbatim logic from the parasite draft) =======================
TP, TN, UN = '#e8575a', '#7b68ee', '#b0b0b0'
TPi, TNi = '#ff5c5c', '#4fd2ff'
MATCH_PX = 15; TOL = 0.1
ROOT = f'{P}/data/parasite_task_results'; FAST = f'{ROOT}/fast-mode/session-20260925-182436'
CXC = f'{P}/Codex runs/parasite_parasite_full_auto_task12_20260928_161707/workspace/output'
PUB = 'published (Fig. 4C)'
nuclei, _ = gt.load_task1(f'{ROOT}/oneframe_results/ExampleData_TunelMeasurement.xlsx')
rois, _, _ = gt.load_task1_coords(f'{P}/data/one_frame', nuclei)
lab = {n['roi']: n for n in nuclei}; rois = sorted(rois, key=lambda r: r['roi'])
hil1 = {}
with open(f'{ROOT}/hil_results/repeat_1/data/c2_tunel_final/C2_TUNEL_Classified_Nuclei.csv') as fh:
    for r in csv.DictReader(fh): hil1[int(r['roi_name'].split()[-1])] = (float(r['corrected_tunel_signal']), r['group'])
RUN1 = {'automated': f'{ROOT}/auto-mode/task_1', 'fast': f'{FAST}/task_1', 'codex': CXC}
R1 = {m: load_run(f'{d}/tunel_nuclei.csv') for m, d in RUN1.items()}
M1 = {m: tf.imread(f'{d}/tunel_nuclei_mask.tif') for m, d in RUN1.items()}
def nearest(runs, r):
    k, (x, y, v, g) = min(runs.items(), key=lambda kv: (kv[1][0] - r['x']) ** 2 + (kv[1][1] - r['y']) ** 2); return v, g
T1 = {'published': [(lab[r['roi']]['corrected'], lab[r['roi']]['group']) for r in rois], 'hil': [hil1[r['roi']] for r in rois]}
T1.update({m: [nearest(R1[m], r) for r in rois] for m in R1})
agree = {m: sum(a[1] == b[1] for a, b in zip(T1[m], T1['published'])) for m in T1}
r1 = {m: np.corrcoef([v for v, _ in T1['published']], [v for v, _ in T1[m]])[0, 1] for m in T1}
ov = ij_roi.read(f'{P}/data/one_frame/TUNEL_Overview.roi'); ox, oy, ow, oh = ov['left'], ov['top'], ov['w'], ov['h']
rawd = f'{ROOT}/auto-mode/task_1/leishmania_tunel_position40_task1/raw_images'
dapi = stretch(tf.imread(f'{rawd}/position40_DAPI_raw16bit_ch3.tif')[oy:oy + oh, ox:ox + ow], 1.0, 99.99) ** 0.6
tun = stretch(tf.imread(f'{rawd}/position40_TUNEL_raw16bit_ch2.tif')[oy:oy + oh, ox:ox + ow], 50.0, 99.97)

LAB2 = '/parin/task2_lab'
g2 = gt.load_task2_coords(LAB2); w2, _ = gt.load_task2(LAB2)
def hil_t2(path, setc, grpc, valc):
    with open(path) as fh: rows = list(csv.DictReader(fh))
    used, out = set(), []
    for r in rows:
        s_, g, v = r[setc], r[grpc].strip(), float(r[valc])
        cand = [i for i, w in enumerate(w2) if w['group'] == g and i not in used and any(x in s_ for x in w['subsets'])] or \
               [i for i, w in enumerate(w2) if w['group'] == g and i not in used]
        i = min(cand, key=lambda i: abs(w2[i]['index'] - v)); used.add(i); out.append((w2[i]['index'], v, g, r))
    return out
HP = f'{ROOT}/hil_results'
H2 = {'run 1': hil_t2(f'{HP}/repeat_1/data/proliferation_index/Proliferation_Index_Unique_Parasites.csv', 'source_sets', 'display_group', 'proliferation_index_value'),
      'run 2': hil_t2(f'{HP}/repeat_2/data/proliferation/position40_proliferation_index_unique_parasites.csv', 'set_id', 'group', 'proliferation_index')}
for g in g2:
    best = min(H2['run 2'], key=lambda t: (float(t[3]['parasite_fullfield_x']) + float(t[3]['parasite_fullfield_w']) / 2 - g['x']) ** 2 +
                                          (float(t[3]['parasite_fullfield_y']) + float(t[3]['parasite_fullfield_h']) / 2 - g['y']) ** 2)
    g['lab_index'] = best[0]
RUN2 = {'automated': f'{ROOT}/auto-mode/task_2', 'fast': f'{FAST}/task_2', 'codex': CXC}
A2 = {}
for m, d in RUN2.items():
    objs = sa.objects_of(bw.auto_table(d, 2)); kept, miss, extra = sa.match_by_distance(g2, objs, MATCH_PX)
    A2[m] = dict(objs=objs, kept=kept, extra=extra,
                 pairs=[(g['lab_index'], r['value'], g['group'], r['group']) for g, r, _ in kept if r['value'] is not None])
t2, _ = dp.collect(ROOT, f'{P}/data/published/jci.insight.169020.sdval.xlsx')
t2 = {k: t2[k] for k in (PUB, 'Agent + human (run 1)', 'Advanced', 'Fast')}
t2['Codex'] = dp.group_values([dict(index=o['value'], group=o['group']) for o in A2['codex']['objs']])
t2 = {k: t2[k] for k in (PUB, 'Codex', 'Fast', 'Advanced', 'Agent + human (run 1)')}
K2 = list(t2); MK = ['published', 'codex', 'fast', 'automated', 'hil']
def mwu(g):
    a, b = g.get('TUNEL+'), g.get('TUNEL-')
    return mannwhitneyu(a, b, alternative='two-sided').pvalue if a and b and len(a) > 1 and len(b) > 1 else None
p2 = [mwu(t2[k]) for k in K2]
S2 = {}
for k, v in H2.items(): S2['hil ' + k] = dict(found=15, extra=0, grp=15, err=[abs(a_ - b_) for a_, b_, _, _ in v])
for m in A2: S2[m] = dict(found=len(A2[m]['kept']), extra=len(A2[m]['extra']), grp=sum(gg == rg for _, _, gg, rg in A2[m]['pairs']),
                          err=[abs(a_ - b_) for a_, b_, _, _ in A2[m]['pairs']])

# ======================= a / e  the two draft pipeline tables (steps + time + verdict) =======================
def pfmt(p): return f'p = {p:.2g}'
def mark(ax, x, y, key, ok):
    ax.text(x, y, '✓' if ok == 'ok' else '✗', ha='center', va='center', fontsize=7, fontweight='bold',
            color=COL[key] if ok == 'ok' else BAD, fontfamily='DejaVu Sans')
def cells(ax, x0, cw, gap, y, rh, cl):
    for i, (tt, kind) in enumerate(cl):
        xx = x0 + i * (cw + gap)
        ax.add_patch(Rectangle((xx, y), cw, rh, fc=FILL[kind], ec='k', lw=0.35))
        ax.text(xx + cw / 2, y + rh / 2, tt, ha='center', va='center', fontsize=FS, linespacing=1.0)
# --- a: CD177 (verbatim content of the CD177 draft, panel a)
def block_title(y, bold, rest):                 # experiment name (bold, 6 pt) + one-line question (5 pt)
    t = fig.text(6 / W, 1 - y / H, bold, fontsize=FL, fontweight='bold', va='top')
    w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi * 25.4
    fig.text((6 + w * 1.08 + 3) / W, 1 - (y + 0.35) / H, rest, fontsize=FS, va='top')   # measured width runs a little short
letter(0, 0, 'a'); block_title(0.3, 'CD177 neutrophil morphology', 'Do CD177+ and CD177− neutrophils differ in cell and nucleus morphology?')
TA_Y, TA_H = 3.6, 36.5; a1 = ax_mm(0, TA_Y, W, TA_H); a1.axis('off'); a1.set_xlim(0, W); a1.set_ylim(TA_H, 0)
steps = ['Cell\nsegmentation', 'Nucleus\nsegmentation', 'QC &\npost-processing', 'Features', 'Statistics']
rows = [('Published', 'published', [('Cellpose (default)', 'manual'), ('ilastik pixel\nclassifier', 'manual'), ('dilation, crop,\nsize filter', 'manual'), ('25 hand-picked', 'manual'), ('cell-level\nMann–Whitney', 'manual')]),
        ('Codex CLI', 'codex', [('Cellpose cyto3 at\nhalf resolution', 'other'), ('erosion/threshold\nwithin cell', 'other'), ('shape QC', 'other'), ('30 features', 'other'), ('donor-level\nWilcoxon, FDR', 'other')]),
        ('Exploratory\nAgentic-J', 'fast', [('Cellpose cyto3', 'auto'), ('Otsu threshold\nwithin cell', 'auto'), ('border, size,\nsolidity filter', 'auto'), ('25 features', 'auto'), ('donor-level\nWilcoxon, FDR', 'auto')]),
        ('Conventional\nAgentic-J', 'automated', [('Cellpose cyto3\non optical density', 'auto'), ('threshold\nwithin cell', 'auto'), ('size/shape QC', 'auto'), ('32 features', 'auto'), ('donor-level\npermutation, FDR', 'auto')]),
        ('Conventional\nAgentic-J\n+ human', 'hil', [('Cellpose-SAM,\nchosen by the user', 'hf'), ('micro-sam fine-tuned on\n11 hand-annotated tiles', 'hf'), ('lobe-to-cell\nassignment', 'auto'), ('11 features', 'auto'), ('image-level\nMann–Whitney, FDR', 'auto')])]
TIME = {'published': (3 * 7 * 24 * 60, 'weeks'), 'hil': (488, '488 min'), 'automated': (67, '67 min'), 'fast': (35, '35 min'), 'codex': (13, '13 min')}
VERD = {'published': ('ref', 'no meaningful difference'), 'hil': ('ok', '0/11 features significant'),
        'automated': ('ok', '0/32 significant, 3/3 runs'), 'fast': ('x', '5/25 significant, 3/3 runs'), 'codex': ('ok', '0/30 features significant')}
x0, cw, gap, rh, y0 = 13, 19.4, 0.5, 5.4, 6.0
STEP_END = x0 + 5 * (cw + gap) - gap
for i, s_ in enumerate(steps): a1.text(x0 + i * (cw + gap) + cw / 2, 3.1, s_, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
TX, VX = STEP_END + 1.2, STEP_END + 33.5
a1.text(TX + 14, 3.1, 'Time to result', ha='center', va='center', fontsize=FS, fontweight='bold')
a1.text(VX + 16, 3.1, 'Scientific verdict', ha='center', va='center', fontsize=FS, fontweight='bold')
for r, (name, key, cl) in enumerate(rows):
    y = y0 + r * (rh + 0.5)
    a1.text(0.3, y + rh / 2, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=0.95)
    cells(a1, x0, cw, gap, y, rh, cl)
    ok, txt = VERD[key]
    if ok != 'ref': mark(a1, VX + 1, y + rh / 2, key, ok)
    a1.text(VX + 3, y + rh / 2, txt, va='center', fontsize=FS, color=COL['published'] if ok == 'ref' else 'k')
ta = ax_mm(TX, TA_Y + y0, 28, 5 * (rh + 0.5) - 0.5)
for r, (name, key, _) in enumerate(rows):
    v, lab_ = TIME[key]; ta.barh(r, v, height=0.55, color=COL[key], lw=0); ta.text(v * 1.3, r, lab_, va='center', fontsize=FS)
ta.set_xscale('log'); ta.set_xlim(5, 4e5); ta.set_ylim(4.5, -0.5); ta.set_yticks([]); ta.minorticks_off()
ta.set_xticks([10, 60, 1440, 10080]); ta.set_xticklabels(['10 min', '1 h', '1 d', '1 wk']); ta.spines['left'].set_visible(False)
# --- shared legend for both tables
YL = TA_Y + TA_H + 3.4; al = ax_mm(0, YL, W, 3); al.axis('off'); al.set_xlim(0, W); al.set_ylim(3, 0)   # below the time-axis labels
xl = 13; rend = fig.canvas.get_renderer()
for k, t in [('manual', 'Manual (analyst)'), ('hf', 'Human feedback (user)'), ('auto', 'Automatic (Agentic-J)'), ('other', 'Automatic (other agent)')]:
    al.add_patch(Rectangle((xl, 0.4), 2.4, 2.2, fc=FILL[k], ec='k', lw=0.3)); tx = al.text(xl + 3.1, 1.5, t, va='center', fontsize=FS)
    xl += 3.1 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0
for ok, t in (('ok', 'published conclusion reproduced'), ('x', 'not reproduced')):
    al.text(xl + 1, 1.5, '✓' if ok == 'ok' else '✗', ha='center', va='center', fontsize=7, fontweight='bold',
            color='#444444' if ok == 'ok' else BAD, fontfamily='DejaVu Sans')
    tx = al.text(xl + 3, 1.5, t, va='center', fontsize=FS); xl += 3 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0

# ======================= b  CD177 gallery (all 4 rows of the CD177 draft) =======================
XL, XR = 0.0, 92.0                          # left column = CD177, right column = Leishmania (as in the scorecard)
yb = YL + 5.0; letter(XL, yb, 'b'); head(XL + 6, yb + 0.6, 'CD177: segmentations of four example cells')
clab = ['raw', 'Published', 'Codex', 'Exploratory\nAgentic-J', 'Conventional\nAgentic-J', 'Conventional\nAgentic-J\n+ human']
ccol = ['k', COL['published'], COL['codex'], COL['fast'], COL['automated'], COL['hil']]
pk = pd.read_csv(f'{CD}/inputs/gallery_picks_main.csv')
def rowlab(stem, ses): d = stem.replace('-', '_').split('_'); return f"{d[0]} {d[1]}\n{ses}"
GX, GXR, tw, tg, rpit = XL + 6.0, XR + 9.5, 12.4, 0.4, 12.9   # b has no row names, e keeps Field / Cell labels
GRID_W = 6 * (tw + tg) - tg                  # image-grid width; the overlay legends live inside it
for ri in range(4):
    p_ = pk.iloc[ri]
    for ci in range(6):
        t = plt.imread(f'/work/inputs/gallery_tiles_main/r{ri}_c{ci}.png')   # re-rendered: black 'missed' outline
        if ri == 0 and ci == 2:                 # Codex: small false nucleus fragment at the cell's left edge (tile px x 209-244, y 439-546)
            ARR = (215, 495)
        ax = ax_mm(GX + ci * (tw + tg), yb + 10.5 + ri * rpit, tw, tw); ax.imshow(t, interpolation='none'); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri == 0: ax.set_title(clab[ci], fontsize=FS, color=ccol[ci], pad=1.5, linespacing=0.95)
        if ri == 0 and ci == 2:
            ax.annotate('', xy=ARR, xytext=(ARR[0] - 150, ARR[1] + 230), arrowprops=dict(arrowstyle='-|>', color='w', lw=0.9, mutation_scale=6))
lx = ax_mm(GX, yb + 10.5 + 4 * rpit + 0.6, GRID_W, 2.6); lx.axis('off'); lx.set_xlim(0, GRID_W); lx.set_ylim(0, 1)
for x, c, ls, t in [(0, '#f2c200', '-', 'cell'), (10, '#27c7f0', '-', 'nucleus'), (23.5, '#FF8C00', '--', 'segmented, not measured'),
                    (51.5, '#000000', ':', 'published cell missed')]:
    lx.plot([x, x + 3.5], [0.5, 0.5], color=c, ls=ls, lw=1); lx.text(x + 4.5, 0.5, t, va='center', fontsize=FS)
lx2 = ax_mm(GX, yb + 10.5 + 4 * rpit + 3.4, GRID_W, 2.6); lx2.axis('off'); lx2.set_xlim(0, GRID_W); lx2.set_ylim(0, 1)
lx2.annotate('', xy=(3.6, 0.5), xytext=(0, 0.5), arrowprops=dict(arrowstyle='-|>', color='#555555', lw=0.9, mutation_scale=6))
lx2.text(4.5, 0.5, 'Codex: small fragment segmented as a nucleus, not a real nucleus', va='center', fontsize=FS)

# ======================= c  CD177 agreement with the published masks =======================
yc = yb; letter(XR, yc, 'c')
ag = pd.read_csv(f'{CD}/inputs/mask_agreement_per_cell.csv')
ag = ag[~ag.method.isin(['hil', 'hil3'])].replace({'method': {'hil4': 'hil'}})
ag['matched'] = (ag.cell_iou >= 0.5) & (ag.measured.astype(str) == 'True')
head(XR + 6, yc + 0.6, f'CD177: agreement with the published masks, {int((ag.method == "hil").sum())} cells')
M = ['codex', 'fast', 'automated', 'hil']
mets = [('matched', 'Cells found & measured,\nIoU ≥ 0.5'), ('cell_iou', 'Cell overlap,\nIoU'), ('nuc_dice', 'Nucleus overlap,\nDice')]
for j, (k, t) in enumerate(mets):
    d = ax_mm(XR + 7 + j * 27.0, yc + 9.5, 23.0, 13.5)
    if k == 'matched':
        for mi, m in enumerate(M): d.bar(mi, ag[ag.method == m].matched.mean(), 0.62, color=COL[m], lw=0)
        d.set_ylim(0, 1.0); d.set_yticks([0, 0.5, 1.0])
    else:
        data = [ag[(ag.method == m) & ag.matched][k].dropna().values for m in M]
        bp = d.boxplot(data, positions=range(4), widths=0.55, whis=(5, 95), showfliers=False, patch_artist=True,
                       medianprops=dict(color='k', lw=0.8), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.5))
        for b, m in zip(bp['boxes'], M): b.set_facecolor(COL[m]); b.set_alpha(0.75)
        d.set_ylim(0.4, 1.0); d.set_yticks([0.4, 0.6, 0.8, 1.0])
    d.set_xlim(-0.6, 3.6); d.set_title(t, fontsize=FS, pad=2, linespacing=0.95); d.set_xticks([])
    d.grid(axis='y', lw=0.25, color='#dddddd'); d.set_axisbelow(True)
    if j == 2: d.set_yticklabels([])
fig.legend(handles=[Patch(fc=COL[m], label=SH[m]) for m in M], frameon=False, fontsize=FS, ncol=2, loc='upper left',
           bbox_to_anchor=((XR + 6) / W, 1 - (yc + 25.0) / H), handlelength=1.0, handletextpad=0.4, columnspacing=1.0, borderaxespad=0)

# ======================= d  CD177 per-cell distributions vs published (verbatim from the CD177 draft, panel d) =======================
ye = yc + 31.5; letter(XR, ye, 'd'); head(XR + 6, ye + 0.6, 'CD177: per-cell distributions vs published data')
PX2 = 0.2515 ** 2
pcs = pd.read_csv(f'{CD}/inputs/published_cells_from_crops.csv'); pcs = pcs[~pcs.is60x]
ref = pd.DataFrame({'cond': pcs.cond, 'cell': pcs.cell_px * PX2, 'nuc': pcs.nuc_px * PX2, 'nc': pcs.nc})
A = pd.read_csv('/auto_cells/cd177_cell_features.csv'); A = A[A.size_analysis_eligible.astype(str) == 'True']
F = pd.read_csv('/fast_cells/cell_features.csv'); F = F[~F.image.astype(str).str.contains('60x')]
Cx = pd.read_csv('/codex_cells/cd177_cell_features.csv'); Cx = Cx[~Cx.magnification.astype(str).str.contains('60')]
Hh = pd.read_csv('/hil_cells/Cell_Morphology_Measurements_V2.csv'); Hh = Hh[Hh.included_primary.astype(str) == 'True']
def std(df, cond, c, n, r): return pd.DataFrame({'cond': np.where(df[cond].astype(str).str.lower().str.contains('pos'), 'pos', 'neg'), 'cell': df[c], 'nuc': df[n], 'nc': df[r]})
DATA = {'published': ref, 'automated': std(A, 'condition', 'cell_area', 'nucleus_area', 'nuclear_fraction_of_cell'),
        'fast': std(F, 'condition', 'cell_area', 'nucleus_area', 'nuclear_fraction'), 'hil': std(Hh, 'CD177_status', 'cell_area_um2', 'nucleus_union_area_um2', 'nuclear_to_cell_area_ratio'),
        'codex': std(Cx, 'condition', 'cell_area', 'nucleus_area', 'nucleus_cell_area_ratio')}
def fmt_p(p):                                  # -> (base, exponent); exponent drawn as its own 5 pt text, raised
    if p >= 0.01: return (f'p = {p:.2f}', None)
    m_, e_ = f'{p:.1e}'.split('e'); return (f'p = {m_}×10', f'{int(e_)}'.replace('-', '−'))
_rend_p = fig.canvas.get_renderer()
def put_p(ax, parts, color):
    """p-value at the axes' top right: base 6 pt, exponent 5 pt raised (mathtext would print it at 3.5 pt)."""
    base, ex = parts
    if ex is None:
        ax.text(1.0, 1.0, base, transform=ax.transAxes, ha='right', va='top', fontsize=6, color=color); return
    e = ax.text(1.0, 1.0, ex, transform=ax.transAxes, ha='right', va='top', fontsize=5, color=color)
    w = e.get_window_extent(_rend_p).width / ax.get_window_extent(_rend_p).width
    ax.text(1.0 - w, 0.84, base, transform=ax.transAxes, ha='right', va='top', fontsize=6, color=color)
MWU = {}
for m in ['automated', 'fast', 'hil', 'codex']:
    for k in ['cell', 'nuc', 'nc']:
        v = DATA[m]; a_ = v[v.cond == 'pos'][k].dropna(); b_ = v[v.cond == 'neg'][k].dropna()
        MWU[(m, k)] = mannwhitneyu(a_, b_, alternative='two-sided').pvalue
PUBP = (('p = 3.5×10', '−52'), ('p = 0.2', None), ('p = 2.4×10', '−11'))      # the paper's own values (as in the CD177 draft)
FEAT = [('cell', 'Cell area (µm²)', (0, 450)), ('nuc', 'Nucleus area (µm²)', (0, 160)), ('nc', 'Nucleus / cell area', (0, 1))]
RM = ['published', 'codex', 'fast', 'automated', 'hil']
RL = {'published': 'Published', 'automated': 'Conventional\nAgentic-J', 'fast': 'Exploratory\nAgentic-J', 'hil': 'Conventional\nAgentic-J\n+ human', 'codex': 'Codex'}
NEG, POS = '#f08080', '#7b68ee'
for ri, m in enumerate(RM):
    for fj, (k, t, lim) in enumerate(FEAT):
        h = ax_mm(XR + 15 + fj * 24.4, ye + 5.5 + ri * 5.5, 22.4, 4.6); v = DATA[m]; bins = np.linspace(*lim, 46)
        for cnd, colr in [('neg', NEG), ('pos', POS)]:
            h.hist(v[v.cond == cnd][k].dropna().clip(*lim), bins=bins, density=True, color=colr, alpha=0.55, lw=0)
        h.set_ylim(0, h.get_ylim()[1] * 2.0)   # headroom: p-value label clear of the bars
        h.set_xlim(lim); h.set_xticks({'cell': [0, 200, 400], 'nuc': [0, 50, 100], 'nc': [0, 0.5, 1]}[k]); h.set_yticks([]); h.spines['left'].set_visible(False)
        txt, tc = (PUBP[fj], COL['published']) if m == 'published' else (fmt_p(MWU[(m, k)]), '#333333')
        put_p(h, txt, tc)
        if ri < len(RM) - 1: h.set_xticklabels([])
        else: h.set_xlabel(t)
        if fj == 0: h.text(-0.05, 0.4, RL[m], transform=h.transAxes, ha='right', va='center', fontsize=FS, fontweight='bold', color=COL[m], linespacing=0.95)
fig.legend(handles=[Patch(fc=NEG, alpha=0.55, label='CD177−'), Patch(fc=POS, alpha=0.55, label='CD177+')], frameon=False, fontsize=FS, ncol=2,
           loc='upper right', bbox_to_anchor=((XR + 88) / W, 1 - (ye + 0.3) / H), handlelength=1.0, columnspacing=0.8, borderaxespad=0)

YP = max(yb + 10.5 + 4 * rpit + 6.4, ye + 5.5 + 5 * 5.5 + 6.6)   # lower block starts under the taller CD177 column
# --- e: parasite (content of the parasite draft, panel a; verdict wording as the Fig. 2 scorecard)
letter(0, YP, 'e'); block_title(YP + 0.3, 'Leishmania infection', 'Task 1: call each of 6 host nuclei TUNEL+, meaning apoptotic, or TUNEL−.   Task 2: compare the proliferation index of 15 parasites between TUNEL+ and TUNEL− hosts.')
YE = YP + 3.6; TE_H = 41.5; a2 = ax_mm(0, YE, W, TE_H); a2.axis('off'); a2.set_xlim(0, W); a2.set_ylim(TE_H, 0)
steps = ['Registration,\nlive → fixed', 'Host / nucleus\nsegmentation', 'TUNEL call', 'Parasite\nlocalisation', 'Measurement\n& statistics']
rows = [('Published', 'published', [('by hand', 'manual'), ('hand-drawn\nROIs', 'manual'), ('by inspection', 'manual'), ('15 ROIs + 3\nbackgrounds\neach, by hand', 'manual'), ('ROI Manager →\nExcel;\nMann–Whitney', 'manual')]),
        ('Codex CLI', 'codex', [('SIFT + RANSAC,\nECC-refined\nto 1.6 px', 'other'), ('T1 DAPI\nthreshold; T2\nCellpose cyto3', 'other'), ('largest-gap\ncutoff', 'other'), ('red-channel\npeaks + water-\nshed: 10 kept', 'other'), ('annulus\nbackground;\nno test', 'other')]),
        ('Exploratory\nAgentic-J', 'fast', [('TurboReg rigid\non brightfield', 'auto'), ('T1 StarDist;\nT2 Cellpose\ncyto3, 7 hosts', 'auto'), ('largest-gap\ncutoff', 'auto'), ('C2 + C3 spot\nevidence:\n22 kept', 'auto'), ('annulus\nbackground;\nno test', 'auto')]),
        ('Conventional\nAgentic-J', 'automated', [('TurboReg /\nStackReg rigid', 'auto'), ('T1 DoG + water-\nshed; T2 Cellpose\ncyto3, 6 hosts', 'auto'), ('2-component\nGMM', 'auto'), ('micro-sam\nvit_b_lm:\n105 → 5 kept', 'auto'), ('annulus\nbackground;\nno test', 'auto')]),
        ('Conventional\nAgentic-J\n+ human', 'hil', [('lab ROIs given,\norigin\ntranslated', 'hf'), ('lab ROIs,\ngiven', 'hf'), ('largest-gap\ncutoff', 'hf'), ('lab ROIs paired\n1:1 with\nbackgrounds', 'hf'), ('native 16-bit;\nMann–Whitney', 'hf')])]
PTIME = {'published': 'not recorded', 'hil': '34 / 26 min', 'automated': '55 / 85 min', 'fast': '22 / 69 min', 'codex': '18 min,\nboth tasks'}
w2n = sum(e <= TOL for e in S2['hil run 1']['err']); sc_ = S2['codex']; pc = p2[K2.index('Codex')]
VER = {'published': [('ref', '2/6 nuclei TUNEL+'), ('ref', f'8 TUNEL+ vs 7 TUNEL−, {pfmt(p2[K2.index(PUB)])}')],
       'hil': [('ok', f'{agree["hil"]}/6 calls, r = {r1["hil"]:.2f}'), ('ok', f'15/15; {w2n}/15 within {TOL}; {pfmt(p2[K2.index("Agent + human (run 1)")])}')],
       'automated': [('ok', f'{agree["automated"]}/6 calls, r = {r1["automated"]:.2f}'), ('x', f'{S2["automated"]["found"]}/15 found, {S2["automated"]["extra"]} false, no TUNEL+')],
       'fast': [('ok', f'{agree["fast"]}/6 calls, r = {r1["fast"]:.2f}, 1 split'), ('x', f'{S2["fast"]["found"]}/15 found, {S2["fast"]["extra"]} false, no host')],
       'codex': [('ok', f'{agree["codex"]}/6 calls, r = {r1["codex"]:.2f}'),
                 ('ok' if sc_['grp'] == sc_['found'] else 'x', f'{sc_["found"]}/15 found, hosts right' + (f'; {pfmt(pc)}' if pc else ''))]}
def t1_text(m):
    """'k/6 calls agree, r = ..' + which lab nuclei differ / were split, so a 5/6 is not read as a failure."""
    diff = [rois[i]['roi'] for i in range(6) if T1[m][i][1] != T1['published'][i][1]]
    split = []
    if m in R1:
        cnt = {}
        for x, y, _v, _g in R1[m].values():
            n = min(rois, key=lambda r: (r['x'] - x) ** 2 + (r['y'] - y) ** 2)['roi']; cnt[n] = cnt.get(n, 0) + 1
        split = sorted(n for n, c in cnt.items() if c > 1)
    notes = []
    if diff: notes.append(('only cell ' if len(diff) == 1 else 'cells ') + ', '.join(map(str, diff)) + (' differs' if len(diff) == 1 else ' differ'))
    if split: notes.append('cell ' + ', '.join(map(str, split)) + ' split in two')
    return f'{6 - len(diff)}/6 calls agree, r = {r1[m]:.2f}' + ('\n' + '; '.join(notes) if notes else '')
for m_ in ('codex', 'fast', 'automated', 'hil'): VER[m_][0] = ('ok', t1_text(m_))
x0, cw, gap, rh, y0 = 14.5, 15.0, 0.5, 6.2, 6.4
STEP_END = x0 + 5 * (cw + gap) - gap
for i, s_ in enumerate(steps): a2.text(x0 + i * (cw + gap) + cw / 2, 3.2, s_, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
xt = STEP_END + 1.0; vx = [xt + 15.0, xt + 15.0 + 33.5]
a2.text(xt + 7, 3.2, 'Time\nT1 / T2', ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
for j, t in enumerate(['Task 1  TUNEL+ / TUNEL− calls\n6 nuclei, fixed frame', 'Task 2  proliferation index\n15 parasites, last live frame']):
    a2.text(vx[j] + 16, 3.2, t, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
for r, (name, key, cl) in enumerate(rows):
    y = y0 + r * (rh + 0.5)
    a2.text(0.3, y + rh / 2, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=0.95)
    cells(a2, x0, cw, gap, y, rh, cl)
    a2.text(xt + 7, y + rh / 2, PTIME[key], ha='center', va='center', fontsize=FS, color='#555555' if key == 'published' else 'k', linespacing=1.0)
    for j, (ok, txt) in enumerate(VER[key]):
        if ok != 'ref': mark(a2, vx[j] + 1, y + rh / 2, key, ok)
        a2.text(vx[j] + 3, y + rh / 2, txt, va='center', fontsize=FS, color=COL['published'] if ok == 'ref' else 'k', linespacing=1.0)

# ======================= f  parasite task-1 gallery (field, cell 3, cell 6), same columns as b =======================
yd = YE + TE_H + 1.5; letter(XL, yd, 'f'); head(XL + 6, yd + 0.6, 'Leishmania task 1: six nuclei, called by each pipeline')
LBL = 5.0                                   # in-image labels (nucleus numbers, scale bars): Nature minimum
def draw_calls(ax, kind):
    if kind in ('published', 'hil'):
        for r, (v, g) in zip(rois, T1[kind]):
            ax.add_patch(Ellipse((r['x'] - ox, r['y'] - oy), r['w'], r['h'], fill=False, lw=0.6, ec=TPi if g == 'TUNEL+' else TNi))
    else:
        for k, (x, y, v, g) in R1[kind].items():
            for c in find_contours((M1[kind] == int(k)).astype(float), 0.5):
                ax.plot(c[:, 1] - ox, c[:, 0] - oy, lw=0.6, color=TPi if g == 'TUNEL+' else TNi)
cols = ['raw', 'published', 'codex', 'fast', 'automated', 'hil']
c3 = [r for r in rois if r['roi'] == 3][0]; c6 = [r for r in rois if r['roi'] == 6][0]
views = [('Field', (0, ow, 0, oh)), ('Cell 3', (c3['x'] - ox - 45, c3['x'] - ox + 45, c3['y'] - oy - 45, c3['y'] - oy + 45)),
         ('Cell 6', (c6['x'] - ox - 45, c6['x'] - ox + 45, c6['y'] - oy - 45, c6['y'] - oy + 45))]
for ri, (rl, (xa, xb, ya, yb_)) in enumerate(views):
    for ci, kind in enumerate(cols):
        ax = ax_mm(XL + 9.5 + ci * (tw + tg), yd + 10.5 + ri * rpit, tw, tw)
        if kind == 'raw': ax.imshow(tun, cmap='magma', interpolation='nearest')
        else: ax.imshow(dapi, cmap='gray', interpolation='nearest'); draw_calls(ax, kind)
        ax.set_xlim(xa, xb); ax.set_ylim(yb_, ya); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri == 0: ax.set_title('raw\nTUNEL' if ci == 0 else clab[ci], fontsize=FS, color=ccol[ci], pad=1.5, linespacing=0.95)
        if ri == 0 and kind != 'raw':
            if kind in ('published', 'hil'): pts = [(r['x'], r['y'], str(r['roi']), 'right') for r in rois]
            else:
                near = {}
                for x, y, _v, _g in R1[kind].values(): near.setdefault(min(rois, key=lambda r: (r['x'] - x) ** 2 + (r['y'] - y) ** 2)['roi'], []).append((x, y))
                pts = []
                for n, v in near.items():
                    v = sorted(v)
                    pts += [(x, y, f'{n}{"ab"[j]}' if len(v) > 1 else str(n), 'left' if (len(v) > 1 and j == 0) else 'right') for j, (x, y) in enumerate(v)]
            for x, y, n, side in pts:
                if side == 'left': ax.text(x - ox + 2, y - oy + 16, n, color='w', fontsize=LBL, ha='center', va='top')
                else: ax.text(x - ox + 20, y - oy + (34 if n == '1' else -6), n, color='w', fontsize=LBL)
        if ci == 0:
            ax.set_ylabel(rl, fontsize=FS, rotation=0, ha='right', va='center', labelpad=1.5)
            if ri < 2:
                L, lab_ = (20, '20 µm') if ri == 0 else (5, '5 µm'); px = L / gt.PIXEL_UM; span = xb - xa
                ax.plot([xa + 0.06 * span, xa + 0.06 * span + px], [yb_ - 0.08 * span] * 2, color='w', lw=0.9)
                ax.text(xa + 0.06 * span, yb_ - 0.12 * span, lab_, color='w', fontsize=LBL, va='bottom')
lx = ax_mm(XL + 9.5, yd + 10.5 + 3 * rpit + 0.6, GRID_W, 2.6); lx.axis('off'); lx.set_xlim(0, GRID_W); lx.set_ylim(0, 1)
for x, c, t in [(0, TPi, 'labeled TUNEL+'), (20, TNi, 'labeled TUNEL−')]:
    lx.plot([x, x + 3.5], [0.5, 0.5], color=c, lw=1); lx.text(x + 4.5, 0.5, t, va='center', fontsize=FS)

# ======================= f  parasite task 1: corrected TUNEL signal per nucleus (from the parasite draft) =======================
yf = yd + 10.5 + 3 * rpit + 5.0; letter(XL, yf, 'g'); head(XL + 6, yf + 0.6, 'Leishmania task 1: corrected TUNEL signal per nucleus')
s = ax_mm(XL + 11, yf + 9.5, 77, 13.0)        # legend sits above the axes, clear of the data
MARK = {'published': 'o', 'hil': 'D', 'automated': 's', 'fast': '^', 'codex': 'v'}
xs = np.arange(6); ML = ['published', 'codex', 'fast', 'automated', 'hil']
for m, off in zip(ML, np.linspace(-0.2, 0.2, 5)):
    v = [a_[0] for a_ in T1[m]]; g = [a_[1] for a_ in T1[m]]
    s.scatter(xs + off, v, marker=MARK[m], s=9, color=COL[m], lw=0, zorder=3, label=SH[m])
    hit = [i for i in range(6) if g[i] == 'TUNEL+']
    s.scatter(xs[hit] + off, np.array(v)[hit], marker=MARK[m], s=34, facecolors='none', edgecolors=TP, lw=0.7, zorder=4)
s.set_xticks(xs); s.set_xticklabels([f'Cell {r["roi"]}' for r in rois]); s.set_xlim(-0.5, 5.5)
s.set_ylabel('Corrected TUNEL\nsignal (a.u.)', linespacing=0.95); s.set_ylim(0, 31); s.set_yticks([0, 10, 20, 30])
s.grid(axis='y', lw=0.25, color='#dddddd'); s.set_axisbelow(True)
h_, l_ = s.get_legend_handles_labels()
h_.append(Line2D([], [], marker='o', ms=4, lw=0, mfc='none', mec=TP, mew=0.7)); l_.append('labeled TUNEL+')
s.legend(h_, l_, frameon=False, fontsize=FS, loc='lower left', bbox_to_anchor=(-0.02, 1.02), ncol=3, handletextpad=0.3,
         borderaxespad=0, labelspacing=0.25, columnspacing=1.0)

# ======================= h  count comparison with the published analysis (= parasite draft c) =======================
yh = yd; letter(XR, yh, 'h'); head(XR + 6, yh + 0.6, 'Leishmania: counts compared with the published analysis')
SCM = [('codex', 'codex', 'Codex CLI'), ('fast', 'fast', 'Exploratory\nAgentic-J'), ('automated', 'automated', 'Conventional\nAgentic-J'),
       ('hil run 1', 'hil', 'Conventional\nAgentic-J\n+ human')]
def t1val(m): return agree['hil' if m.startswith('hil') else m]
def within(m): return sum(e <= TOL for e in S2[m]['err'])
SROWS = [('Task 1', 'TUNEL calls agree', lambda m: (t1val(m), 6)),
         ('Task 2', 'Lab parasites found', lambda m: (S2[m]['found'], 15)),
         ('', 'Host group correct', lambda m: (S2[m]['grp'], S2[m]['found'])),
         ('', f'Index within {TOL} of lab', lambda m: (within(m), len(S2[m]['err']))),
         ('', 'False detections', lambda m: S2[m]['extra'])]
sc = ax_mm(XR, yh + 3.5, 88, 26); sc.axis('off'); sc.set_xlim(0, 88); sc.set_ylim(26, 0)
LX, CX0, CW, CG, RY0, RH, RG = 0.3, 34.5, 12.8, 0.6, 6.4, 3.1, 0.4
for j, (m, cm, lab_) in enumerate(SCM):
    sc.text(CX0 + j * (CW + CG) + CW / 2, 3.0, lab_, ha='center', va='center', fontsize=FS, color=COL[cm], linespacing=0.95, fontweight='bold')
y = RY0
for i, (grp, lab_, fn) in enumerate(SROWS):
    if grp == 'Task 2': y += 1.0
    if grp: sc.text(LX, y + RH / 2, grp, va='center', fontsize=FS, fontweight='bold')
    sc.text(LX + 8.5, y + RH / 2, lab_, va='center', fontsize=FS, color='#555555' if lab_ == 'False detections' else 'k')
    for j, (m, cm, _) in enumerate(SCM):
        x = CX0 + j * (CW + CG); v = fn(m)
        if lab_ == 'False detections':
            sc.add_patch(Rectangle((x, y), CW, RH, fc='#f2f2f2', ec='none'))
            sc.text(x + CW / 2, y + RH / 2, str(v), ha='center', va='center', fontsize=FS, color='#333333'); continue
        k, n = v
        if n == 0:
            sc.add_patch(Rectangle((x, y), CW, RH, fc='w', ec='#bbbbbb', lw=0.3, ls='--'))
            sc.text(x + CW / 2, y + RH / 2, '–', ha='center', va='center', fontsize=FS, color='#888888'); continue
        fr = k / n; sc.add_patch(Rectangle((x, y), CW, RH, fc=COL[cm], alpha=0.12 + 0.88 * fr, ec='none'))
        sc.text(x + CW / 2, y + RH / 2, f'{k} / {n}', ha='center', va='center', fontsize=FS, color='w' if fr > 0.55 else 'k',
                fontweight='bold' if fr > 0.55 else 'normal')
    y += RH + RG

# ======================= i  parasite task 2 = published Fig. 4C (+ delta median, from the parasite draft) =======================
yg = yh + 30.5; letter(XR, yg, 'i'); head(XR + 6, yg + 0.6, 'Leishmania task 2: proliferation index by host TUNEL status')
SRC = ['Published', 'Codex', 'Exploratory', 'Conventional', 'Conventional\n+ human']
def stars(p): return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
rng = np.random.default_rng(0)
EH = H - 1.5 - 15.3 - 9.0 - yg; d = ax_mm(XR + 11, yg + 9.0, 77, EH); lim = (0, 10.6)   # plot fills the column to the page bottom; lim = (0, 10.6)   # legend sits above, plot fills the column
for i, k in enumerate(K2):
    g = t2[k]
    if i % 2 == 0: d.axvspan(i - 0.5, i + 0.5, color='#f1f1f1', lw=0, zorder=0)
    for grp, dx, c in [('TUNEL+', -0.19, TP), ('TUNEL-', 0.19, TN), ('Unassigned', 0, UN)]:
        v = g.get(grp)
        if not v: continue
        n_out = sum(x_ < lim[0] for x_ in v)
        d.scatter(i + dx + rng.uniform(-0.05, 0.05, len(v)), np.clip(v, *lim), s=6, color=c, lw=0, alpha=0.85, zorder=2)
        d.plot([i + dx - 0.14, i + dx + 0.14], [st.median(v)] * 2, color='k', lw=0.8, zorder=4, solid_capstyle='butt')
        if n_out: d.text(i + dx, 0.15, f'↓{n_out}', ha='center', va='bottom', fontsize=LBL, color=c)
    if p2[i] is not None: d.text(i, lim[1], stars(p2[i]), ha='center', va='top', fontsize=FS)
d.set_xticks(range(len(K2))); d.set_xticklabels(SRC, fontsize=FS, rotation=40, ha='right', rotation_mode='anchor', linespacing=0.95)
for tl, m in zip(d.get_xticklabels(), MK): tl.set_color(COL[m])
d.tick_params(axis='x', length=0); d.set_ylim(*lim); d.set_yticks([0, 5, 10]); d.set_xlim(-0.5, len(K2) - 0.5)
d.grid(axis='y', lw=0.25, color='#dddddd'); d.set_axisbelow(True); d.set_ylabel('Proliferation index,\n10 − R/G', linespacing=0.95)
d.legend(handles=[Line2D([], [], marker='o', lw=0, ms=2.5, color=TP, label='TUNEL+ host'), Line2D([], [], marker='o', lw=0, ms=2.5, color=TN, label='TUNEL− host'),
                  Line2D([], [], marker='o', lw=0, ms=2.5, color=UN, label='no host'), Line2D([], [], color='k', lw=0.8, label='median')],
         frameon=False, fontsize=FS, loc='lower left', bbox_to_anchor=(-0.02, 1.03), ncol=4, handlelength=0.9, handletextpad=0.3, borderaxespad=0, columnspacing=1.2)
def dmed(g, B=4000):
    a_, b_ = np.asarray(g.get('TUNEL+', [])), np.asarray(g.get('TUNEL-', []))
    if len(a_) < 2 or len(b_) < 2: return None
    r = np.random.default_rng(1); bs = [np.median(r.choice(a_, len(a_))) - np.median(r.choice(b_, len(b_))) for _ in range(B)]
    return np.median(a_) - np.median(b_), *np.percentile(bs, [2.5, 97.5])
DM = [dmed(t2[k]) for k in K2]
import matplotlib.transforms as mtr
ty = mtr.blended_transform_factory(d.transData, fig.transFigure); yrow = 1 - (yg + 9.0 + EH + 11.6) / H; yci = 1 - (yg + 9.0 + EH + 13.8) / H
d.text(-0.5, yrow, 'Δ median  ', ha='right', va='center', fontsize=FS, fontweight='bold', transform=ty)
d.text(-0.5, yci, '95% CI  ', ha='right', va='center', fontsize=FS, color='#555555', transform=ty)
for i, r in enumerate(DM):
    npos = len(t2[K2[i]].get('TUNEL+', []))
    d.text(i, yrow, (f'{r[0]:.2f}' if r is not None else ('no\nTUNEL+' if npos == 0 else f'only {npos}\nTUNEL+')), ha='center', va='center',
           fontsize=FS, color=COL[MK[i]], style='normal' if r is not None else 'italic', linespacing=0.9, transform=ty)
    if r is not None: d.text(i, yci, f'[{r[1]:.1f}, {r[2]:.1f}]', ha='center', va='center', fontsize=FS, color='#555555', transform=ty)

for ext in ['pdf', 'svg', 'png']: fig.savefig(f'/work/Fig2_full.{ext}', dpi=600 if ext == 'png' else None)   # 600 dpi: combination artwork
print('CD177 block: left ends', round(yb + 10.5 + 4 * rpit + 6.0, 1), 'right ends', round(ye + 5.5 + 5 * 5.5 + 6, 1), '| Leishmania: left ends', round(yf + 9.5 + 13.0 + 3.5, 1), 'right ends', round(yg + 9.0 + EH + 15.3, 1), '| i plot height', round(EH, 1), '| page', H)
# what was plotted, so a rebuild can be checked against the drafts
print('task1 calls agree', agree, ' r', {m: round(v, 3) for m, v in r1.items()})
for m, v in S2.items(): print('task2', f'{m:10}', {k: v[k] for k in ('found', 'extra', 'grp')})
print('task2 p', {k: (None if p is None else float(f'{p:.3g}')) for k, p in zip(K2, p2)})
