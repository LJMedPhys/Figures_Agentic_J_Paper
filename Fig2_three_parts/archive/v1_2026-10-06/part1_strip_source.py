# ======================= part 1: brief plugin cases, a-c (no YOFF: drawn with ax_abs) =======================
from scipy import ndimage
VARIANT = os.environ.get('FIG2_VARIANT', 'orientationj')       # 'orientationj' (placeholder) or 'coloc2'
TW, TG = 58.0, 3.0                                              # tile width, gap between tiles (mm)
TH_Y, TH_S = 5.2, 21.0                                          # thumbnail top and size
PL_X, PL_W, PL_Y, PL_H = 31.5, 25.5, 5.6, 17.5                  # result plot, relative to the tile
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
    ar = ax_abs(x + TH_S + 1.2, TH_Y + TH_S / 2 - 1.5, 5.5, 3); ar.axis('off'); ar.set_xlim(0, 1); ar.set_ylim(0, 1)
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
    tb.text(0.03, 0.97, '250 mM NaCl', transform=tb.transAxes, color='w', fontsize=FS, va='top')
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
for v_ in np.unique(lab0)[1:]: rgb[lab0 == v_] = np.array(cmap_l(int(v_) % 20)[:3]) * 0.55   # cells dimmed so the tracks read on top
tc = thumb(xc, rgb)
spots = list(csv.DictReader(open('/work/inputs/tracking/mouse_cell_tracking_trackmate/mouse_cell_spots.csv')))[3:]
trk = {}
for r_ in spots:
    if r_['TRACK_ID'] not in ('', 'None'):
        trk.setdefault(r_['TRACK_ID'], []).append((int(float(r_['FRAME'])), float(r_['POSITION_X']), float(r_['POSITION_Y'])))
for t_ in trk.values():
    t_.sort(); tc.plot([p[1] for p in t_], [p[2] for p in t_], color='w', lw=0.35, alpha=0.9)
tc.set_xlim(0, lab0.shape[1]); tc.set_ylim(lab0.shape[0], 0)
tc.text(0.03, 0.97, f'{len(trk)} tracks, 101 frames', transform=tc.transAxes, color='w', fontsize=FS, va='top')
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
