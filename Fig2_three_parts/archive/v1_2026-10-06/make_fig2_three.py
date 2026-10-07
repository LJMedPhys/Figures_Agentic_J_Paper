"""Fig. 2, three-part version (v1, 2026-10-06), following the feedback: (1) brief visuals of plugin cases, then
(2) CD177 and (3) Leishmania. 180 mm wide, 231 mm deep (legend on a separate page), text >= 5 pt, PNG 600 dpi.

  Part 1, plugin cases, one small tile each: input image -> result -> one-line takeaway
    a  BoneJ, trabecular bone morphometry, Smith et al. 2023      (inputs: ../BoneJ_reproduction/inputs, inputs/bonej)
    b  OrientationJ placeholder, or with FIG2_VARIANT=coloc2: Coloc 2 colocalization, Simmons et al. 2022
                                                                  (inputs: ../Coloc2_reproduction/inputs)
    c  TrackMate, mouse stem cell tracking                        (inputs: inputs/tracking, unpacked from
                                                                   ../tracking_mouse_stem_cells.zip)
  Parts 2 and 3: the content of ../Fig2_validation/make_fig2.py v10, moved down by YOFF and relettered:
    d  scorecard: CD177 frame on the left = part 2, Leishmania frame on the right = part 3
    e-g  CD177: gallery, cells found and measured, per-cell distributions
    h-j  Leishmania: gallery, TUNEL signal per nucleus, proliferation index

Generated once from ../Fig2_validation/make_fig2.py v10; the two scripts do not track each other afterwards.
run.sh builds both variants. Mounts: as ../Fig2_validation/run.sh, plus /fv = ../Fig2_validation/inputs,
/bonej = ../BoneJ_reproduction/inputs, /coloc = ../Coloc2_reproduction/inputs (all read-only), /work = this folder."""
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
MM = 1 / 25.4; P1H = 41.0; YOFF = P1H; W, H = 180, 190 + P1H   # part 1 strip on top
fig = plt.figure(figsize=(W * MM, H * MM))
def ax_mm(x, y, w, h, **k): return fig.add_axes([x / W, 1 - (YOFF + y + h) / H, w / W, h / H], **k)
def ax_abs(x, y, w, h, **k): return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **k)   # part 1, no offset
def letter(x, y, s): fig.text(x / W, 1 - (YOFF + y) / H, s, fontsize=8, fontweight='bold', va='top')
def head(x, y, s): fig.text(x / W, 1 - (YOFF + y) / H, s, fontsize=FL, va='top')
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

# ======================= part 1: brief plugin cases, a-c (no YOFF: drawn with ax_abs) =======================
from scipy import ndimage
VARIANT = os.environ.get('FIG2_VARIANT', 'orientationj')       # 'orientationj' (placeholder) or 'coloc2'
TW, TG = 58.0, 3.0                                              # tile width, gap between tiles (mm)
TH_Y, TH_S = 5.2, 19.5                                          # thumbnail top and size
PL_X, PL_W, PL_Y, PL_H = 34.5, 23.0, 5.6, 17.5                  # result plot, relative to the tile (clear of the arrow)
MSG_Y = 30.4                                                    # takeaway line(s)
_r1 = fig.canvas.get_renderer()
GREEN, BLUE, PURPLE = '#1a9641', '#0571b0', '#7b3294'

def tile(i, lab, plugin, task):
    x = i * (TW + TG)
    fig.text(x / W, 1, lab, fontsize=8, fontweight='bold', va='top')
    t = fig.text((x + 5) / W, 1 - 0.35 / H, plugin, fontsize=FL, fontweight='bold', va='top')
    w = t.get_window_extent(_r1).width / fig.dpi * 25.4
    fig.text((x + 5 + w + 1.5) / W, 1 - 0.55 / H, task, fontsize=FS, va='top')
    return x
def arrow(x):                                                   # input -> result
    ar = ax_abs(x + TH_S + 0.8, TH_Y + TH_S / 2 - 1.5, 4.2, 3); ar.axis('off'); ar.set_xlim(0, 1); ar.set_ylim(0, 1)
    ar.annotate('', xy=(1, 0.5), xytext=(0, 0.5), arrowprops=dict(arrowstyle='-|>', color='#555555', lw=0.8, mutation_scale=6))
def message(x, txt):
    fig.text(x / W, 1 - MSG_Y / H, txt, fontsize=FS, va='top', linespacing=1.15)
def thumb(x, img, cmap=None, aspect_box=1.0):
    h = TH_S * aspect_box
    ax = ax_abs(x, TH_Y + (TH_S - h) / 2, TH_S, h); ax.imshow(img, cmap=cmap, interpolation='antialiased'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_linewidth(0.3)
    return ax

# ---------- a  BoneJ
xa = tile(0, 'a', 'BoneJ', 'trabecular bone morphometry')
f_ = np.load('/work/inputs/bonej/FMNH_218302_tball_L03_mask.npz'); sh = tuple(f_['shape'])
vol = np.unpackbits(f_['packed'])[:int(np.prod(sh))].reshape(sh)
proj = vol.mean(1).T                                            # side projection through the stack, slices left to right
ta = thumb(xa, proj, cmap='gray', aspect_box=proj.shape[0] / proj.shape[1])
ta.text(0.03, 0.05, 'L3 bone mask', transform=ta.transAxes, color='w', fontsize=FS, va='bottom')
arrow(xa)
ref = {f"{r['genus']} {r['species']}": r for r in csv.DictReader(open('/bonej/reference_table2.csv'))}
summ = {r['species']: r for r in csv.DictReader(open('/bonej/species_summary.csv'))}
BJ = [('BV/TV', 'bone_volume_fraction_mean', 'bv_tv'), ('Tb.Th', 'mean_trabecular_thickness_mm_mean', 'tb_th_mm'),
      ('Tb.Sp', 'mean_trabecular_spacing_mm_mean', 'tb_sp_mm'), ('Conn.D', 'connectivity_density_per_mm3_mean', 'conn_d_per_mm3')]
dev = {lab_: [100 * (float(summ[sp][c]) - float(ref[sp][rc])) / float(ref[sp][rc]) for sp in ref] for lab_, c, rc in BJ}
pa = ax_abs(xa + PL_X, PL_Y, PL_W, PL_H)
pa.axvspan(-1, 1, color='#eeeeee', lw=0, zorder=0); pa.axvline(0, color='#999999', lw=0.5, zorder=1)
rng1 = np.random.default_rng(5)
for k, (lab_, _, _) in enumerate(BJ):
    v = np.array(dev[lab_]); pa.scatter(v, k + rng1.uniform(-0.18, 0.18, len(v)), s=3.5, color=GREEN, lw=0, zorder=3)
pa.set_yticks(range(4)); pa.set_yticklabels([b[0] for b in BJ]); pa.set_ylim(3.6, -0.6)
pa.set_xlim(-2, 12); pa.set_xticks([0, 5, 10]); pa.set_xlabel('Agentic-J − published, %')
pa.tick_params(axis='y', length=0)
ok3 = max(abs(d) for lab_ in ('BV/TV', 'Tb.Th', 'Conn.D') for d in dev[lab_])
sp_lo, sp_hi = min(dev['Tb.Sp']), max(dev['Tb.Sp'])
message(xa, f'66 vertebrae, 11 species: BV/TV, Tb.Th and Conn.D within\n{ok3:.1f}% of published; Tb.Sp +{sp_lo:.0f} to +{sp_hi:.0f}%, from an unreported setting')

# ---------- b  OrientationJ placeholder, or Coloc 2
if VARIANT == 'coloc2':
    xb = tile(1, 'b', 'Coloc 2', 'colocalization under osmotic stress')
    raw = tf.imread('/coloc/examples/nacl_field1_raw.tif').astype(float); mk = tf.imread('/coloc/examples/nacl_field1_mask.tif') > 0
    S_ = 660; cov = ndimage.uniform_filter(mk.astype(float), S_, mode='constant'); h2 = S_ // 2
    cy, cx = np.unravel_index(np.argmax(cov[h2:-h2, h2:-h2]), cov[h2:-h2, h2:-h2].shape); sl = (slice(cy, cy + S_), slice(cx, cx + S_))
    nrm = lambda a: np.clip((a - np.percentile(a, 0.5)) / (np.percentile(a, 99.8) - np.percentile(a, 0.5)), 0, 1)
    red, grn = nrm(raw[..., 0][sl]), nrm(raw[..., 1][sl])
    tb = thumb(xb, np.stack([red, grn, red], -1))
    tb.text(0.03, 0.97, '250 mM NaCl', transform=tb.transAxes, color='w', fontsize=FS, va='top',
            bbox=dict(boxstyle='square,pad=0.2', fc='k', ec='none', alpha=0.65))
    tb.text(0.03, 0.05, 'Su(Hw)', transform=tb.transAxes, color='#ff6ad5', fontsize=FS, va='bottom', fontweight='bold')
    tb.text(0.52, 0.05, 'γH2Av', transform=tb.transAxes, color='#7dff7d', fontsize=FS, va='bottom', fontweight='bold')
    arrow(xb)
    def ld(name, col):
        rows = list(csv.DictReader(open(f'/coloc/{name}')))
        return (np.array([float(r[col]) for r in rows if r['condition'].startswith('iso')]),
                np.array([float(r[col]) for r in rows if r['condition'].startswith('hyper')]))
    paper, agentc = ld('paper_fig5_pcc.csv', 'pearson_r_above_threshold'), ld('v2_interactive_coloc2.csv', 'pearson_r_above_costes')
    pb = ax_abs(xb + PL_X, PL_Y, PL_W, PL_H)
    rng2 = np.random.default_rng(3)
    for ci_, x0 in enumerate((0.0, 1.2)):
        for dx, vals, col in ((0.0, paper[ci_], PURPLE), (0.42, agentc[ci_], BLUE)):
            pb.scatter(x0 + dx + rng2.uniform(-0.08, 0.08, len(vals)), vals, s=3.5, color=col, lw=0, zorder=3)
            pb.plot([x0 + dx - 0.15, x0 + dx + 0.15], [vals.mean()] * 2, color='k', lw=0.8, zorder=4)
    pb.axhline(0, color='#999999', lw=0.4, ls='--', zorder=0)
    pb.set_xticks([0.21, 1.41]); pb.set_xticklabels(['isotonic', 'NaCl']); pb.tick_params(axis='x', length=0)
    pb.set_ylim(-0.45, 0.85); pb.set_yticks([0, 0.4, 0.8]); pb.set_ylabel("Pearson's R", labelpad=1)
    pb.text(0.02, 1.0, 'published', transform=pb.transAxes, color=PURPLE, fontsize=FS, va='bottom')
    pb.text(0.52, 1.0, 'Agentic-J', transform=pb.transAxes, color=BLUE, fontsize=FS, va='bottom')
    d_pub, d_agt = paper[1].mean() - paper[0].mean(), agentc[1].mean() - agentc[0].mean()
    message(xb, f'Osmotic stress raises nuclear γH2Av–Su(Hw) correlation:\nΔ mean R +{d_agt:.2f}, published +{d_pub:.2f}, 17 image fields')
else:
    xb = tile(1, 'b', 'OrientationJ', 'orientation analysis')
    from matplotlib.patches import FancyBboxPatch as _FB
    ph = ax_abs(xb, TH_Y, TW, MSG_Y - TH_Y + 5.5); ph.axis('off'); ph.set_xlim(0, 1); ph.set_ylim(0, 1)
    ph.add_patch(_FB((0.01, 0.02), 0.98, 0.96, boxstyle='round,pad=0,rounding_size=0.04', fc='#f7f7f7', ec='#aaaaaa', lw=0.6, ls='--'))
    ph.text(0.5, 0.56, 'Placeholder', ha='center', va='center', fontsize=FL, fontweight='bold', color='#777777')
    ph.text(0.5, 0.40, 'OrientationJ task: input image, result, one-line takeaway', ha='center', va='center', fontsize=FS, color='#777777')

# ---------- c  TrackMate
xc = tile(2, 'c', 'TrackMate', 'mouse stem cell tracking')
lab0 = tf.imread('/work/inputs/tracking/mouse_cell_tracking/man_track1275.tif')
cmap_l = mpl.colormaps['tab20']
rgb = np.zeros(lab0.shape + (3,))
rgb[lab0 > 0] = 0.42                                            # cells of the first frame in grey; tracks in colour on top
tc = thumb(xc, rgb)
spots = list(csv.DictReader(open('/work/inputs/tracking/mouse_cell_tracking_trackmate/mouse_cell_spots.csv')))[3:]
trk = {}
for r_ in spots:
    if r_['TRACK_ID'] not in ('', 'None'):
        trk.setdefault(r_['TRACK_ID'], []).append((int(float(r_['FRAME'])), float(r_['POSITION_X']), float(r_['POSITION_Y'])))
for k_, t_ in enumerate(sorted(trk.values(), key=len, reverse=True)):
    t_.sort(); tc.plot([p[1] for p in t_], [p[2] for p in t_], color=cmap_l(k_ % 20), lw=0.4, alpha=0.95)
tc.set_xlim(0, lab0.shape[1]); tc.set_ylim(lab0.shape[0], 0)
tc.text(0.03, 0.97, f'{len(trk)} tracks', transform=tc.transAxes, color='w', fontsize=FS, va='top',
        bbox=dict(boxstyle='square,pad=0.2', fc='k', ec='none', alpha=0.65))
arrow(xc)
iv = [float(r_['mean_velocity_um_per_min']) for r_ in csv.DictReader(open('/work/inputs/tracking/plot_extra_modifications/interval_mean_velocities.csv'))]
pc = ax_abs(xc + PL_X, PL_Y, PL_W, PL_H)
rng3 = np.random.default_rng(7)
pc.scatter(rng3.uniform(-0.28, 0.28, len(iv)), iv, s=3.0, color='#444444', lw=0, zorder=3)
pc.plot([-0.38, 0.38], [np.mean(iv)] * 2, color='k', lw=0.9, zorder=4)
pc.text(0.42, np.mean(iv), f'{np.mean(iv):.1f}', va='center', fontsize=FS)
pc.set_xlim(-0.6, 0.9); pc.set_xticks([]); pc.set_ylim(0, 3.5); pc.set_yticks([0, 1, 2, 3])
pc.set_ylabel('Mean speed per frame\ninterval, µm/min', linespacing=0.95, labelpad=1)
pc.spines['bottom'].set_visible(False)
n_spots = sum(len(t_) for t_ in trk.values())
message(xc, f'{len(trk)} tracks, {n_spots} spots over 101 frames, 5 min apart;\nmean speed {np.mean(iv):.1f} µm/min across {len(iv)} frame intervals')

# separator between part 1 and parts 2-3
sep = ax_abs(0, P1H - 3.2, W, 0.5); sep.axis('off'); sep.set_xlim(0, W); sep.set_ylim(0, 1); sep.axhline(0.5, color='#bbbbbb', lw=0.5)
print('part 1:', VARIANT, '| BoneJ max |dev| 3 metrics', round(ok3, 2), 'Tb.Sp', round(sp_lo, 1), round(sp_hi, 1),
      '| tracks', len(trk), 'spots', n_spots, 'mean interval speed', round(float(np.mean(iv)), 3))

# ======================= a  scorecard =======================
letter(0, 0, 'd'); a = ax_mm(0, 0, W, 44); a.axis('off'); a.set_xlim(0, W); a.set_ylim(44, 0)
ROWS = [('Published', 'published', ['manual']), ('Codex CLI', 'codex', ['other']),
        ('Exploratory\nAgentic-J', 'fast', ['auto']), ('Conventional\nAgentic-J', 'automated', ['auto']),
        ('Conventional\nAgentic-J\n+ human', 'hil', ['hf', 'auto'])]
# CD177: values as hard-coded in the CD177 draft (v7)
CD_TIME = {'published': 'weeks', 'codex': '13 min', 'fast': '35 min', 'automated': '67 min', 'hil': '488 min'}
CD_VER = {'published': ('ref', 'no meaningful difference'), 'codex': ('plain', '0/30 features significant'),
          'fast': ('plain', '5/25 significant, 3/3 runs'), 'automated': ('plain', '0/32 significant, 3/3 runs'),
          'hil': ('plain', '0/11 tested significant\n34 features extracted')}
# Biologist feedback on the segmentations (Lukas, 2026-10-06)
CD_FB = {'published': ('ref', 'reference segmentation'),
         'codex': ('plain', 'many cells missed;\nfaulty cell and nucleus segmentation'),
         'fast': ('plain', 'faulty nucleus segmentation;\nfalse-positive cells'),
         'automated': ('plain', 'cells segmented but discarded;\nnuclei missed or partly segmented'),
         'hil': ('plain', 'clean cell segmentation;\nnuclei rounded and over-segmented')}
# parasite: recomputed above, same wording as the parasite draft (shortened to fit)
PA_TIME = {'published': 'not recorded', 'codex': '18 min,\nboth tasks', 'fast': '22 / 69 min', 'automated': '55 / 85 min', 'hil': '34 / 26 min'}
def pfmt(p): return f'p = {p:.2g}'
w2n = sum(e <= TOL for e in S2['hil run 1']['err']); sc_ = S2['codex']; pc = p2[K2.index('Codex')]
def t1_text(m):
    """The run's own calls, listed as TUNEL+ / TUNEL- nuclei like the published entry (as in Fig2_full)."""
    pos = [r['roi'] for r, (_v, g) in zip(rois, T1[m]) if g == 'TUNEL+']
    neg = [r['roi'] for r, (_v, g) in zip(rois, T1[m]) if g != 'TUNEL+']
    split = []
    if m in R1:
        cnt = {}
        for x, y, _v, _g in R1[m].values():
            n = min(rois, key=lambda r: (r['x'] - x) ** 2 + (r['y'] - y) ** 2)['roi']; cnt[n] = cnt.get(n, 0) + 1
        split = sorted(n for n, c in cnt.items() if c > 1)
    def cl(v): return ('nucleus ' if len(v) == 1 else 'nuclei ') + ', '.join(f'{n}, split in two' if n in split else str(n) for n in v)
    return f'TUNEL+: {cl(pos)}\nTUNEL−: {cl(neg)}'
PA_T1 = {'published': ('ref', 'TUNEL+: nuclei 3, 6\nTUNEL−: nuclei 1, 2, 4, 5')}
PA_T1.update({m_: ('plain', t1_text(m_)) for m_ in ('codex', 'fast', 'automated', 'hil')})
PA_T2 = {'published': ('ref', f'higher proliferation in TUNEL+ hosts\n8 TUNEL+ vs 7 TUNEL−, {pfmt(p2[K2.index(PUB)])}'),
         'codex': ('plain', 'higher proliferation in TUNEL+ hosts\n' + (f'{pfmt(pc)}; ' if pc else '') + f'{sc_["found"]} of 15 parasites'),
         'fast': ('plain', f'not testable: no parasite found\n{S2["fast"]["extra"]} false detections, no host'),
         'automated': ('plain', f'not testable: no TUNEL+ host found\n{S2["automated"]["found"]} of 15 parasites, {S2["automated"]["extra"]} false'),
         'hil': ('plain', f'higher proliferation in TUNEL+ hosts\n{pfmt(p2[K2.index("Agent + human (run 1)")])}; all 15 parasites, {w2n} within {TOL}')}

XCHIP = 16.3                                  # analysis chip; frames start right after it
GRP = [(24.5, 95.5, 'CD177 neutrophil morphology'), (97.5, 179.6, 'Leishmania infection: TUNEL and proliferation')]
XCDT, HC_CD, HC_FB = 29.5, 46.5, 77.5         # CD177: time, verdict, biologist feedback (column centres)
XPAT, HC_T1, HC_T2 = 103.5, 124.0, 159.8      # Leishmania: time, task 1, task 2
y0, rh, rp = 9.4, 5.4, 6.0
GBOT = y0 + 5 * rp + 0.1                      # frames end just below the last row
for x0_, x1_, t in GRP:                       # one outlined frame + shaded header band per experiment
    a.add_patch(FancyBboxPatch((x0_, 0.25), x1_ - x0_, GBOT - 0.25, boxstyle='round,pad=0,rounding_size=1.2',
                               fc='none', ec='#8c8c8c', lw=0.6, zorder=3))
    a.add_patch(FancyBboxPatch((x0_, 0.25), x1_ - x0_, 3.0, boxstyle='round,pad=0,rounding_size=1.2',
                               fc='#e6e6e6', ec='none', zorder=1))
    a.add_patch(Rectangle((x0_, 2.0), x1_ - x0_, 1.25, fc='#e6e6e6', ec='none', zorder=1))   # square off the band's bottom
    a.text((x0_ + x1_) / 2, 1.75, t, ha='center', va='center', fontsize=FL, fontweight='bold', zorder=4)
SUB = [(XCHIP + 2.4, 'Analysis'), (XCDT, 'Time'), (HC_CD, 'Scientific verdict'), (HC_FB, 'Biologist feedback'),
       (XPAT, 'Time\nT1 / T2'), (HC_T1, 'Task 1  TUNEL+ and TUNEL− nuclei\n6 nuclei, fixed frame'), (HC_T2, 'Task 2  proliferation in TUNEL+ hosts\n15 parasites, last live frame')]
for x, t in SUB: a.text(x, 6.3, t, ha='center', va='center', fontsize=FS, fontweight='bold', linespacing=1.0)
def verdict(x, y, key, ok, txt):              # v10: entries without tick marks or crosses
    a.text(x, y, txt, va='center', fontsize=FS, color=COL['published'] if ok == 'ref' else 'k', linespacing=1.0)
_rend = fig.canvas.get_renderer()
def _w(txt):                                    # rendered width of an entry in mm
    t = a.text(0, 0, txt, fontsize=FS); w = t.get_window_extent(_rend).width / fig.dpi * 25.4; t.remove(); return w
def block_x(centre, entries):                   # left x so that the widest entry is centred on the heading
    return centre - max(_w(t) for _, t in entries.values()) / 2
X_CD, X_FB, X_T1, X_T2 = block_x(HC_CD, CD_VER), block_x(HC_FB, CD_FB), block_x(HC_T1, PA_T1), block_x(HC_T2, PA_T2)
for r, (name, key, kinds) in enumerate(ROWS):
    y = y0 + r * rp; yc = y + rh / 2
    if r % 2 == 1:                           # stripes break at the gaps, so the two experiments read as separate blocks
        for sx0, sx1 in [(0, GRP[0][0] - 1.0)] + [(g0 + 0.3, g1 - 0.3) for g0, g1, _ in GRP]:
            a.add_patch(Rectangle((sx0, y - 0.3), sx1 - sx0, rh + 0.6, fc='#f4f4f4', ec='none', zorder=0))
    a.text(0.3, yc, name, va='center', fontsize=FS, fontweight='bold', color=COL[key], linespacing=0.95)
    cw = 4.8 / len(kinds)                       # analysis chip: who did the work (split for human + automatic)
    for j, k in enumerate(kinds): a.add_patch(Rectangle((XCHIP + j * cw, yc - 1.6), cw, 3.2, fc=FILL[k], ec='k', lw=0.3))
    a.text(XCDT, yc, CD_TIME[key], ha='center', va='center', fontsize=FS, color='#555555' if key == 'published' else 'k')
    verdict(X_CD, yc, key, *CD_VER[key]); verdict(X_FB, yc, key, *CD_FB[key])
    a.text(XPAT, yc, PA_TIME[key], ha='center', va='center', fontsize=FS, color='#555555' if key == 'published' else 'k', linespacing=1.0)
    verdict(X_T1, yc, key, *PA_T1[key]); verdict(X_T2, yc, key, *PA_T2[key])
xl = XCHIP; rend = fig.canvas.get_renderer(); yl = y0 + 5 * rp + 2.7
for k, t in [('manual', 'Manual (analyst)'), ('hf', 'Human feedback (user)'), ('auto', 'Automatic (Agentic-J)'), ('other', 'Automatic (other agent)')]:
    a.add_patch(Rectangle((xl, yl - 1.1), 2.4, 2.2, fc=FILL[k], ec='k', lw=0.3)); tx = a.text(xl + 3.1, yl, t, va='center', fontsize=FS)
    xl += 3.1 + tx.get_window_extent(rend).width / fig.dpi * 25.4 + 3.0

# ======================= b  CD177 gallery (all 4 rows of the CD177 draft) =======================
XL, XR = 0.0, 92.0                          # left column = CD177, right column = Leishmania (as in the scorecard)
yb = 47; letter(XL, yb, 'e'); head(XL + 6, yb + 0.6, 'CD177: segmentations of four example cells')
clab = ['raw', 'Published', 'Codex', 'Exploratory\nAgentic-J', 'Conventional\nAgentic-J', 'Conventional\nAgentic-J\n+ human']
ccol = ['k', COL['published'], COL['codex'], COL['fast'], COL['automated'], COL['hil']]
pk = pd.read_csv(f'{CD}/inputs/gallery_picks_main.csv')
def rowlab(stem, ses): d = stem.replace('-', '_').split('_'); return f"{d[0]} {d[1]}\n{ses}"
GX, GXR, tw, tg, rpit = XL + 6.0, XR + 9.5, 12.4, 0.4, 12.9   # b has no row names, e keeps Field / Cell labels
GRID_W = 6 * (tw + tg) - tg                  # image-grid width; the overlay legends live inside it
for ri in range(4):
    p_ = pk.iloc[ri]
    for ci in range(6):
        t = plt.imread(f'/fv/gallery_tiles_main/r{ri}_c{ci}.png')   # re-rendered: black 'missed' outline (as Fig2_full)
        ax = ax_mm(GX + ci * (tw + tg), yb + 10.5 + ri * rpit, tw, tw); ax.imshow(t, interpolation='none'); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.3)
        if ri == 0: ax.set_title(clab[ci], fontsize=FS, color=ccol[ci], pad=1.5, linespacing=0.95)
        if ri == 0 and ci == 2:                 # Codex: small false nucleus fragment at the cell's left edge (tile px x 209-244, y 439-546)
            ARR = (215, 495)
            ax.annotate('', xy=ARR, xytext=(ARR[0] - 150, ARR[1] + 230), arrowprops=dict(arrowstyle='-|>', color='#e8000b', lw=1.1, mutation_scale=7))
lx = ax_mm(GX, yb + 10.5 + 4 * rpit + 0.6, GRID_W, 2.6); lx.axis('off'); lx.set_xlim(0, GRID_W); lx.set_ylim(0, 1)
for x, c, ls, t in [(0, '#f2c200', '-', 'cell'), (10, '#27c7f0', '-', 'nucleus'), (23.5, '#FF8C00', '--', 'segmented, not measured'),
                    (51.5, '#000000', ':', 'published cell missed')]:
    lx.plot([x, x + 3.5], [0.5, 0.5], color=c, ls=ls, lw=1); lx.text(x + 4.5, 0.5, t, va='center', fontsize=FS)
lx2 = ax_mm(GX, yb + 10.5 + 4 * rpit + 3.4, GRID_W, 2.6); lx2.axis('off'); lx2.set_xlim(0, GRID_W); lx2.set_ylim(0, 1)
lx2.annotate('', xy=(3.6, 0.5), xytext=(0, 0.5), arrowprops=dict(arrowstyle='-|>', color='#e8000b', lw=1.1, mutation_scale=7))
lx2.text(4.5, 0.5, 'Codex: small fragment segmented as a nucleus, not a real nucleus', va='center', fontsize=FS)

# ======================= c  CD177 cells found & measured (IoU / Dice boxplots: supplementary/Supp_CD177_mask_overlap) =======================
yc = yb + 10.5 + 4 * rpit + 8.3; letter(XL, yc, 'f')
ag = pd.read_csv(f'{CD}/inputs/mask_agreement_per_cell.csv')
ag = ag[~ag.method.isin(['hil', 'hil3'])].replace({'method': {'hil4': 'hil'}})
ag['matched'] = (ag.cell_iou >= 0.5) & (ag.measured.astype(str) == 'True')
head(XL + 6, yc + 0.6, f'CD177: published cells found and measured, n = {int((ag.method == "hil").sum())}')
M = ['codex', 'fast', 'automated', 'hil']
d = ax_mm(XL + 15, yc + 6.5, 30, 15.5)
for mi, m in enumerate(M):
    fr_ = ag[ag.method == m].matched.mean()
    d.bar(mi, fr_, 0.62, color=COL[m], lw=0); d.text(mi, fr_ + 0.03, f'{fr_:.2f}', ha='center', va='bottom', fontsize=FS)
d.set_ylim(0, 1.12); d.set_yticks([0, 0.5, 1.0]); d.set_xlim(-0.6, 3.6); d.set_xticks([])
d.set_ylabel('Fraction found & measured,\nIoU ≥ 0.5', linespacing=0.95)
d.grid(axis='y', lw=0.25, color='#dddddd'); d.set_axisbelow(True)
fig.legend(handles=[Patch(fc=COL[m], label=SH[m]) for m in M], frameon=False, fontsize=FS, ncol=1, loc='center left',
           bbox_to_anchor=((XL + 49) / W, 1 - (YOFF + yc + 14.2) / H), handlelength=1.0, handletextpad=0.4, labelspacing=0.45, borderaxespad=0)

# ======================= d  CD177 per-cell distributions vs published (verbatim from the CD177 draft, panel d) =======================
ye = yc + 25.5; letter(XL, ye, 'g'); head(XL + 6, ye + 0.6, 'CD177: per-cell distributions vs published data')
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
        h = ax_mm(XL + 15 + fj * 24.4, ye + 5.5 + ri * 6.9, 22.4, 5.85); v = DATA[m]; bins = np.linspace(*lim, 46)   # v9: taller rows, c freed the height
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
           loc='upper right', bbox_to_anchor=((XL + 88) / W, 1 - (YOFF + ye + 0.3) / H), handlelength=1.0, columnspacing=0.8, borderaxespad=0)

# ======================= e  parasite task-1 gallery (field, cell 3, cell 6), same columns as b =======================
yd = yb; letter(XR, yd, 'h'); head(XR + 6, yd + 0.6, 'Leishmania task 1: six nuclei, called by each pipeline')
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
views = [('Field', (0, ow, 0, oh)), ('Nucleus 3', (c3['x'] - ox - 45, c3['x'] - ox + 45, c3['y'] - oy - 45, c3['y'] - oy + 45)),
         ('Nucleus 6', (c6['x'] - ox - 45, c6['x'] - ox + 45, c6['y'] - oy - 45, c6['y'] - oy + 45))]
for ri, (rl, (xa, xb, ya, yb_)) in enumerate(views):
    for ci, kind in enumerate(cols):
        ax = ax_mm(GXR + ci * (tw + tg), yd + 10.5 + ri * rpit, tw, tw)
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
lx = ax_mm(GXR, yd + 10.5 + 3 * rpit + 0.6, GRID_W, 2.6); lx.axis('off'); lx.set_xlim(0, GRID_W); lx.set_ylim(0, 1)
for x, c, t in [(0, TPi, 'labeled TUNEL+'), (20, TNi, 'labeled TUNEL−')]:
    lx.plot([x, x + 3.5], [0.5, 0.5], color=c, lw=1); lx.text(x + 4.5, 0.5, t, va='center', fontsize=FS)

# ======================= f  parasite task 1: corrected TUNEL signal per nucleus (from the parasite draft) =======================
yf = yd + 10.5 + 3 * rpit + 5.5; letter(XR, yf, 'i'); head(XR + 6, yf + 0.6, 'Leishmania task 1: corrected TUNEL signal per nucleus')
s = ax_mm(XR + 11, yf + 9.5, 77, 20.5)        # legend sits above the axes, clear of the data
MARK = {'published': 'o', 'hil': 'D', 'automated': 's', 'fast': '^', 'codex': 'v'}
xs = np.arange(6); ML = ['published', 'codex', 'fast', 'automated', 'hil']
for m, off in zip(ML, np.linspace(-0.2, 0.2, 5)):
    v = [a_[0] for a_ in T1[m]]; g = [a_[1] for a_ in T1[m]]
    s.scatter(xs + off, v, marker=MARK[m], s=9, color=COL[m], lw=0, zorder=3, label=SH[m])
    hit = [i for i in range(6) if g[i] == 'TUNEL+']
    s.scatter(xs[hit] + off, np.array(v)[hit], marker=MARK[m], s=34, facecolors='none', edgecolors=TP, lw=0.7, zorder=4)
s.set_xticks(xs); s.set_xticklabels([f'Nucleus {r["roi"]}' for r in rois]); s.set_xlim(-0.5, 5.5)
s.set_ylabel('Corrected TUNEL\nsignal (a.u.)', linespacing=0.95); s.set_ylim(0, 31); s.set_yticks([0, 10, 20, 30])
s.grid(axis='y', lw=0.25, color='#dddddd'); s.set_axisbelow(True)
h_, l_ = s.get_legend_handles_labels()
h_.append(Line2D([], [], marker='o', ms=4, lw=0, mfc='none', mec=TP, mew=0.7)); l_.append('labeled TUNEL+')
s.legend(h_, l_, frameon=False, fontsize=FS, loc='lower left', bbox_to_anchor=(-0.02, 1.02), ncol=3, handletextpad=0.3,
         borderaxespad=0, labelspacing=0.25, columnspacing=1.0)

# ======================= g  parasite task 2 = published Fig. 4C (+ delta median, from the parasite draft) =======================
yg = yf + 34; letter(XR, yg, 'j'); head(XR + 6, yg + 0.6, 'Leishmania task 2: proliferation index by host TUNEL status')
SRC = ['Published', 'Codex', 'Exploratory', 'Conventional', 'Conventional\n+ human']
def stars(p): return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
rng = np.random.default_rng(0)
EH = 28.5; d = ax_mm(XR + 11, yg + 9.0, 77, EH); lim = (0, 10.6)   # legend sits above, plot fills the column
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
d.grid(axis='y', lw=0.25, color='#dddddd'); d.set_axisbelow(True); d.set_ylabel('Proliferation index\n(10 − R/G)', linespacing=0.95)
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
ty = mtr.blended_transform_factory(d.transData, fig.transFigure); yrow = 1 - (YOFF + yg + 9.0 + EH + 11.6) / H; yci = 1 - (YOFF + yg + 9.0 + EH + 13.8) / H
d.text(-0.5, yrow, 'Δ median  ', ha='right', va='center', fontsize=FS, fontweight='bold', transform=ty)
d.text(-0.5, yci, '95% CI  ', ha='right', va='center', fontsize=FS, color='#555555', transform=ty)
for i, r in enumerate(DM):
    npos = len(t2[K2[i]].get('TUNEL+', []))
    d.text(i, yrow, (f'{r[0]:.2f}' if r is not None else ('no\nTUNEL+' if npos == 0 else f'only {npos}\nTUNEL+')), ha='center', va='center',
           fontsize=FS, color=COL[MK[i]], style='normal' if r is not None else 'italic', linespacing=0.9, transform=ty)
    if r is not None: d.text(i, yci, f'[{r[1]:.1f}, {r[2]:.1f}]', ha='center', va='center', fontsize=FS, color='#555555', transform=ty)

for ext in ['pdf', 'svg', 'png']: fig.savefig(f'/work/Fig2_three_parts_{VARIANT}.{ext}', dpi=600 if ext == 'png' else None)
# what was plotted, so a rebuild can be checked against the drafts
print('verdict block x (mm): CD177', round(X_CD, 1), ' T1', round(X_T1, 1), ' T2', round(X_T2, 1))
print('left column ends', round(ye + 5.5 + 5 * 6.9 + 5.5, 1), '| page', H)
print('task1 calls agree', agree, ' r', {m: round(v, 3) for m, v in r1.items()})
for m, v in S2.items(): print('task2', f'{m:10}', {k: v[k] for k in ('found', 'extra', 'grp')})
print('task2 p', {k: (None if p is None else float(f'{p:.3g}')) for k, p in zip(K2, p2)})
