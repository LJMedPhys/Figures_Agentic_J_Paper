"""Fig. 2 (draft v1, 2026-10-01): validation against two published studies, fusing the CD177 and the
parasite draft figures into one display item (Brief Communication: 2 display items in total).

  a  scorecard: the five pipelines x outcomes of both studies (replaces the two step-by-step tables,
     which move to Extended Data)
  b  CD177 gallery, 2 example cells (one Bottom, one Top acquisition batch) x raw + 5 pipelines
  c  CD177 agreement with the published masks (cells found, cell IoU, nucleus Dice)
  d  parasite task-1 gallery (field + cell 3) x raw TUNEL + 5 pipelines, same columns as b
  e  parasite task 2 = published Fig. 4C (proliferation index by host TUNEL status)

Moved out of the main figure (Extended Data): both step tables, CD177 per-cell distributions, parasite
corrected TUNEL signal per nucleus, the parasite agreement table, remaining gallery rows.

Code and numbers are taken from ../CD177_draft/make_figure.py (v7) and ../Parasite_draft/make_figure.py;
those drafts are not modified. CD177 time/verdict values are the ones hard-coded in the CD177 draft;
parasite values are recomputed here from the run outputs exactly as in the parasite draft.
Style: Nature sizes (min 5 pt, headings 6 pt, panel letters 8 pt bold), 180 mm wide.
Mounts (see run.sh): /cd177 = CD177_draft (ro), /p = parasite bundle (ro), /parin = parasite inputs (ro),
/work = this folder."""
import os, sys, csv, statistics as st
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, Ellipse
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
MM = 1 / 25.4; W, H = 180, 133
fig = plt.figure(figsize=(W * MM, H * MM))
def ax_mm(x, y, w, h, **k): return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **k)
def letter(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight='bold', va='top')
def head(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=FL, va='top')
COL = {'published': '#7b3294', 'hil': '#0571b0', 'automated': '#1a9641', 'fast': '#a6611a', 'codex': '#636363'}
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

# ======================= a  scorecard =======================
letter(0, 0, 'a'); a = ax_mm(0, 0, W, 44); a.axis('off'); a.set_xlim(0, W); a.set_ylim(44, 0)
ROWS = [('Published', 'published', ['manual']), ('Codex CLI', 'codex', ['other']),
        ('Exploratory\nAgentic-J', 'fast', ['auto']), ('Conventional\nAgentic-J', 'automated', ['auto']),
        ('Conventional\nAgentic-J\n+ human', 'hil', ['hf', 'auto'])]
# CD177: values as hard-coded in the CD177 draft (v7)
CD_TIME = {'published': 'weeks', 'codex': '13 min', 'fast': '35 min', 'automated': '67 min', 'hil': '488 min'}
CD_VER = {'published': ('ref', 'no meaningful difference'), 'codex': ('ok', '0/30 features significant'),
          'fast': ('x', '5/25 significant (3/3 runs)'), 'automated': ('ok', '0/32 significant (3/3 runs)'),
          'hil': ('ok', '0/11 features significant')}
# parasite: recomputed above, same wording as the parasite draft (shortened to fit)
PA_TIME = {'published': 'not recorded', 'codex': '18 min (both)', 'fast': '22 / 69 min', 'automated': '55 / 85 min', 'hil': '32 min / n.r.'}
def pfmt(p): return f'p = {p:.2g}'
w2n = sum(e <= TOL for e in S2['hil run 1']['err']); sc_ = S2['codex']; pc = p2[K2.index('Codex')]
PA_T1 = {'published': ('ref', '2/6 nuclei TUNEL+'),
         'codex': ('ok', f'{agree["codex"]}/6 calls, r = {r1["codex"]:.2f}'),
         'fast': ('ok', f'{agree["fast"]}/6 calls, r = {r1["fast"]:.2f}, 1 split'),
         'automated': ('ok', f'{agree["automated"]}/6 calls, r = {r1["automated"]:.2f}'),
         'hil': ('ok', f'{agree["hil"]}/6 calls, r = {r1["hil"]:.2f}')}
PA_T2 = {'published': ('ref', f'8 TUNEL+ vs 7 TUNEL−, {pfmt(p2[K2.index(PUB)])}'),
         'codex': ('ok' if sc_['grp'] == sc_['found'] else 'x', f'{sc_["found"]}/15 found, hosts right' + (f'; {pfmt(pc)}' if pc else '')),
         'fast': ('x', f'{S2["fast"]["found"]}/15 found (+{S2["fast"]["extra"]}), no host'),
         'automated': ('x', f'{S2["automated"]["found"]}/15 found (+{S2["automated"]["extra"]}), no TUNEL+'),
         'hil': ('ok', f'15/15; {w2n}/15 within {TOL}; {pfmt(p2[K2.index("Agent + human (run 1)")])}')}

XCHIP, XCDT, XCDV, XPAT, XPT1, XPT2 = 27.5, 43.5, 52.0, 99.0, 107.5, 141.0
GRP = [(36.0, 87.0, 'CD177 neutrophil morphology'), (90.0, 180.0, 'Leishmania infection: TUNEL and proliferation')]
for x0_, x1_, t in GRP:
    a.text((x0_ + x1_) / 2, 1.6, t, ha='center', va='center', fontsize=FL, fontweight='bold')
    a.plot([x0_ + 0.5, x1_ - 0.5], [3.2, 3.2], color='k', lw=0.4)
SUB = [(XCHIP + 2.4, 'Analysis'), (XCDT, 'Time'), (XCDV + 17, 'Scientific verdict'), (XPAT, 'Time\nT1 / T2'),
       (XPT1 + 15, 'Task 1  TUNEL calls\n(6 nuclei)'), (XPT2 + 19, 'Task 2  published Fig. 4C\n(15 parasites)')]
for x, t in SUB: a.text(x, 6.3, t, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
y0, rh, rp = 9.4, 5.4, 6.0
def verdict(x, y, key, ok, txt):
    if ok != 'ref':
        a.text(x + 1, y, '✓' if ok == 'ok' else '✗', ha='center', va='center', fontsize=7, fontweight='bold',
               color=COL[key] if ok == 'ok' else BAD, fontfamily='DejaVu Sans')
    a.text(x + 3, y, txt, va='center', fontsize=FS, color=COL['published'] if ok == 'ref' else 'k')
for r, (name, key, kinds) in enumerate(ROWS):
    y = y0 + r * rp; yc = y + rh / 2
    if r % 2 == 1: a.add_patch(Rectangle((0, y - 0.3), W, rh + 0.6, fc='#f6f6f6', ec='none', zorder=0))
    a.text(0.3, yc, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=0.95)
    cw = 4.8 / len(kinds)                       # analysis chip: who did the work (split for human + automatic)
    for j, k in enumerate(kinds): a.add_patch(Rectangle((XCHIP + j * cw, yc - 1.6), cw, 3.2, fc=FILL[k], ec='k', lw=0.3))
    a.text(XCDT, yc, CD_TIME[key], ha='center', va='center', fontsize=FS, color='#555555' if key == 'published' else 'k')
    verdict(XCDV, yc, key, *CD_VER[key])
    a.text(XPAT, yc, PA_TIME[key], ha='center', va='center', fontsize=FS, color='#555555' if key == 'published' else 'k')
    verdict(XPT1, yc, key, *PA_T1[key]); verdict(XPT2, yc, key, *PA_T2[key])
xl = XCHIP; rend = fig.canvas.get_renderer(); yl = y0 + 5 * rp + 1.6
for k, t in [('manual', 'Manual (analyst)'), ('hf', 'Human feedback (user)'), ('auto', 'Automatic (Agentic-J)'), ('other', 'Automatic (other agent)')]:
    a.add_patch(Rectangle((xl, yl - 1.1), 2.4, 2.2, fc=FILL[k], ec='k', lw=0.3)); tx = a.text(xl + 3.1, yl, t, va='center', fontsize=FS)
    xl += 3.1 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0
a.text(xl + 1, yl, '✓', ha='center', va='center', fontsize=7, fontweight='bold', color='#444444', fontfamily='DejaVu Sans')
tx = a.text(xl + 3, yl, 'published conclusion reproduced', va='center', fontsize=FS)
xl += 3 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0
a.text(xl + 1, yl, '✗', ha='center', va='center', fontsize=7, fontweight='bold', color=BAD, fontfamily='DejaVu Sans')
a.text(xl + 3, yl, 'not reproduced', va='center', fontsize=FS)

# ======================= b  CD177 gallery (rows 0 and 2 of the CD177 draft) =======================
yb = 47; letter(0, yb, 'b'); head(6, yb + 0.6, 'CD177: two example cells (Bottom and Top batch)')
clab = ['raw', 'Published', 'Codex', 'Exploratory\nAgentic-J', 'Conventional\nAgentic-J', 'Conventional\nAgentic-J + human']
ccol = ['k', COL['published'], COL['codex'], COL['fast'], COL['automated'], COL['hil']]
pk = pd.read_csv(f'{CD}/inputs/gallery_picks_main.csv')
def rowlab(stem, ses): d = stem.replace('-', '_').split('_'); return f"{d[0]} {d[1]}\n{ses}"
GX, tw, tg = 12.5, 13.0, 0.4
for ri, src in enumerate([0, 2]):
    p_ = pk.iloc[src]
    for ci in range(6):
        t = plt.imread(f'{CD}/inputs/gallery_tiles_main/r{src}_c{ci}.png')
        ax = ax_mm(GX + ci * (tw + tg), yb + 9.0 + ri * (tw + 0.5), tw, tw); ax.imshow(t, interpolation='none'); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri == 0: ax.set_title(clab[ci], fontsize=FS, color=ccol[ci], pad=1.5, linespacing=0.95)
        if ci == 0: ax.set_ylabel(rowlab(p_.stem, p_.session), fontsize=FS, rotation=0, ha='right', va='center', labelpad=1.5)
lx = ax_mm(GX, yb + 36.3, 80, 2.6); lx.axis('off'); lx.set_xlim(0, 80); lx.set_ylim(0, 1)
for x, c, ls, t in [(0, '#f2c200', '-', 'cell'), (9, '#27c7f0', '-', 'nucleus'), (22, '#FF8C00', '--', 'segmented, not measured'),
                    (52, '#FF1744', ':', 'published cell missed')]:
    lx.plot([x, x + 3.5], [0.5, 0.5], color=c, ls=ls, lw=1); lx.text(x + 4.5, 0.5, t, va='center', fontsize=FS)

# ======================= c  CD177 agreement with the published masks =======================
xc = 97; letter(xc, yb, 'c')
ag = pd.read_csv(f'{CD}/inputs/mask_agreement_per_cell.csv')
ag = ag[~ag.method.isin(['hil', 'hil3'])].replace({'method': {'hil4': 'hil'}})
ag['matched'] = (ag.cell_iou >= 0.5) & (ag.measured.astype(str) == 'True')
head(xc + 6, yb + 0.6, f'CD177: agreement with the published masks ({int((ag.method == "hil").sum())} cells)')
M = ['codex', 'fast', 'automated', 'hil']
mets = [('matched', 'Cells found & measured\n(IoU ≥ 0.5)'), ('cell_iou', 'Cell overlap\n(IoU)'), ('nuc_dice', 'Nucleus overlap\n(Dice)')]
for j, (k, t) in enumerate(mets):
    d = ax_mm(xc + 6 + j * 26.5, yb + 10, 22.5, 15)
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
           bbox_to_anchor=((xc + 5) / W, 1 - (yb + 27.5) / H), handlelength=1.0, handletextpad=0.4, columnspacing=1.0, borderaxespad=0)

# ======================= d  parasite task-1 gallery (field + cell 3), same columns as b =======================
yd = 88; letter(0, yd, 'd'); head(6, yd + 0.6, 'Leishmania task 1: the same six nuclei, called by each pipeline')
def draw_calls(ax, kind):
    if kind in ('published', 'hil'):
        for r, (v, g) in zip(rois, T1[kind]):
            ax.add_patch(Ellipse((r['x'] - ox, r['y'] - oy), r['w'], r['h'], fill=False, lw=0.6, ec=TPi if g == 'TUNEL+' else TNi))
    else:
        for k, (x, y, v, g) in R1[kind].items():
            for c in find_contours((M1[kind] == int(k)).astype(float), 0.5):
                ax.plot(c[:, 1] - ox, c[:, 0] - oy, lw=0.6, color=TPi if g == 'TUNEL+' else TNi)
cols = ['raw', 'published', 'codex', 'fast', 'automated', 'hil']
c3 = [r for r in rois if r['roi'] == 3][0]
views = [('Field', (0, ow, 0, oh)), ('Cell 3', (c3['x'] - ox - 45, c3['x'] - ox + 45, c3['y'] - oy - 45, c3['y'] - oy + 45))]
for ri, (rl, (xa, xb, ya, yb_)) in enumerate(views):
    for ci, kind in enumerate(cols):
        ax = ax_mm(GX + ci * (tw + tg), yd + 6.0 + ri * (tw + 0.5), tw, tw)
        if kind == 'raw': ax.imshow(tun, cmap='magma', interpolation='nearest')
        else: ax.imshow(dapi, cmap='gray', interpolation='nearest'); draw_calls(ax, kind)
        ax.set_xlim(xa, xb); ax.set_ylim(yb_, ya); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri == 0 and ci == 0: ax.set_title('raw TUNEL', fontsize=FS, pad=1.5)
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
                if side == 'left': ax.text(x - ox + 2, y - oy + 14, n, color='w', fontsize=4.2, ha='center', va='top')
                else: ax.text(x - ox + 20, y - oy + (30 if n == '1' else -6), n, color='w', fontsize=4.2)
        if ci == 0:
            ax.set_ylabel(rl, fontsize=FS, rotation=0, ha='right', va='center', labelpad=1.5)
            L, lab_ = (20, '20 µm') if ri == 0 else (5, '5 µm'); px = L / gt.PIXEL_UM; span = xb - xa
            ax.plot([xa + 0.06 * span, xa + 0.06 * span + px], [yb_ - 0.08 * span] * 2, color='w', lw=0.9)
            ax.text(xa + 0.06 * span, yb_ - 0.12 * span, lab_, color='w', fontsize=4.2, va='bottom')
lx = ax_mm(GX, yd + 33.3, 80, 2.6); lx.axis('off'); lx.set_xlim(0, 80); lx.set_ylim(0, 1)
for x, c, t in [(0, TPi, 'labeled TUNEL+'), (20, TNi, 'labeled TUNEL−')]:
    lx.plot([x, x + 3.5], [0.5, 0.5], color=c, lw=1); lx.text(x + 4.5, 0.5, t, va='center', fontsize=FS)

# ======================= e  parasite task 2 = published Fig. 4C =======================
letter(xc, yd, 'e'); head(xc + 6, yd + 0.6, 'Leishmania task 2: proliferation index (published Fig. 4C)')
SRC = ['Published', 'Codex', 'Exploratory', 'Conventional', 'Conventional\n+ human']
def stars(p): return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
rng = np.random.default_rng(0)
d = ax_mm(xc + 10, yd + 6, 50, 27); lim = (0, 10.6)
for i, k in enumerate(K2):
    g = t2[k]
    if i % 2 == 0: d.axvspan(i - 0.5, i + 0.5, color='#f1f1f1', lw=0, zorder=0)
    for grp, dx, c in [('TUNEL+', -0.19, TP), ('TUNEL-', 0.19, TN), ('Unassigned', 0, UN)]:
        v = g.get(grp)
        if not v: continue
        n_out = sum(x_ < lim[0] for x_ in v)
        d.scatter(i + dx + rng.uniform(-0.05, 0.05, len(v)), np.clip(v, *lim), s=6, color=c, lw=0, alpha=0.85, zorder=2)
        d.plot([i + dx - 0.14, i + dx + 0.14], [st.median(v)] * 2, color='k', lw=0.8, zorder=4, solid_capstyle='butt')
        if n_out: d.text(i + dx, 0.15, f'↓{n_out}', ha='center', va='bottom', fontsize=4.2, color=c)
    if p2[i] is not None: d.text(i, lim[1], stars(p2[i]), ha='center', va='top', fontsize=FS)
d.set_xticks(range(len(K2))); d.set_xticklabels(SRC, fontsize=FS, rotation=40, ha='right', rotation_mode='anchor', linespacing=0.95)
for tl, m in zip(d.get_xticklabels(), MK): tl.set_color(COL[m])
d.tick_params(axis='x', length=0); d.set_ylim(*lim); d.set_yticks([0, 5, 10]); d.set_xlim(-0.5, len(K2) - 0.5)
d.grid(axis='y', lw=0.25, color='#dddddd'); d.set_axisbelow(True); d.set_ylabel('Proliferation index\n(10 − R/G)', linespacing=0.95)
d.legend(handles=[Line2D([], [], marker='o', lw=0, ms=2.5, color=TP, label='TUNEL+ host'), Line2D([], [], marker='o', lw=0, ms=2.5, color=TN, label='TUNEL− host'),
                  Line2D([], [], marker='o', lw=0, ms=2.5, color=UN, label='no host'), Line2D([], [], color='k', lw=0.8, label='median')],
         frameon=False, fontsize=FS, loc='upper left', bbox_to_anchor=(1.02, 1.0), handlelength=0.9, handletextpad=0.3, borderaxespad=0, labelspacing=0.35)

for ext in ['pdf', 'svg', 'png']: fig.savefig(f'/work/Fig2_validation.{ext}', dpi=450 if ext == 'png' else None)
# what was plotted, so a rebuild can be checked against the drafts
print('task1 calls agree', agree, ' r', {m: round(v, 3) for m, v in r1.items()})
for m, v in S2.items(): print('task2', f'{m:10}', {k: v[k] for k in ('found', 'extra', 'grp')})
print('task2 p', {k: (None if p is None else float(f'{p:.3g}')) for k, p in zip(K2, p2)})
