"""Coloc 2 reproduction figure, prompt v2 only: Agentic-J vs Simmons et al. 2022 (PLoS Genet), Fig. 5A-C.

Style as Fig. 2: Nature sizes (min 5 pt), 180 mm wide; published = purple, Agentic-J interactive = blue.
Inputs (inputs/): paper_fig5_pcc.csv, v2_interactive_coloc2.csv, examples/ (raw RGB fields + the agent's DAPI masks for
the two fields closest to their group means). Pixel size 0.0645 um (ORCA-ER 6.45 um pixels / 100x; the authors' macro
uses the same value). Outputs: Fig_Coloc2_v2.{pdf,svg,png} and figure_numbers_v2.txt.
"""
import csv, itertools, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import tifffile
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__)); IN = os.path.join(HERE, 'inputs')
FS, FL = 5, 6
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': FS,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'xtick.labelsize': FS, 'ytick.labelsize': FS, 'axes.labelsize': FS, 'axes.titlesize': FS, 'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5})
COL = {'published': '#7b3294', 'hil': '#0571b0'}
GREY = '#bbbbbb'
PX_UM = 0.0645
MM = 1 / 25.4; W, H = 180, 70
fig = plt.figure(figsize=(W * MM, H * MM))
def ax_mm(x, y, w, h, **k): return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **k)
def letter(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight='bold', va='top')
def title(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=FL, va='top')

def load(name, col):
    rows = list(csv.DictReader(open(f'{IN}/{name}')))
    return (np.array([float(r[col]) for r in rows if r['condition'].startswith('iso')]),
            np.array([float(r[col]) for r in rows if r['condition'].startswith('hyper')]))
paper = load('paper_fig5_pcc.csv', 'pearson_r_above_threshold')
agent = load('v2_interactive_coloc2.csv', 'pearson_r_above_costes')

def perm_p(x, y):
    v = np.concatenate([x, y]); obs = abs(y.mean() - x.mean()); n = tot = 0
    for c in itertools.combinations(range(len(v)), len(x)):
        m = np.zeros(len(v), bool); m[list(c)] = True
        n += abs(v[~m].mean() - v[m].mean()) >= obs - 1e-12; tot += 1
    return n / tot
def boot_ci(x, y, n=10000, seed=20261006):
    rng = np.random.default_rng(seed)
    d = rng.choice(y, (n, len(y))).mean(1) - rng.choice(x, (n, len(x))).mean(1)
    return np.percentile(d, [2.5, 97.5])
def fmt_p(p):
    m, e = f'{p:.1e}'.split('e'); return f'{m}×10$^{{{int(e)}}}$'

# ── a: example fields: Su(Hw)/gammaH2Av merge and, separately, DAPI with the agent's nuclear mask ──
letter(0, 0, 'a'); title(5, 0.4, 'Example fields')
EX = [('isotonic_field6', 'Isotonic', agent[0][5]), ('nacl_field1', '250 mM NaCl', agent[1][1])]
raws = [tifffile.imread(f'{IN}/examples/{k}_raw.tif').astype(float) for k, *_ in EX]
masks = [tifffile.imread(f'{IN}/examples/{k}_mask.tif') > 0 for k, *_ in EX]
# one display range per channel, shared by both fields so their intensities stay comparable
lo = [np.percentile(np.concatenate([r[..., c].ravel() for r in raws]), 0.5) for c in (0, 1, 2)]
hi = [np.percentile(np.concatenate([r[..., c].ravel() for r in raws]), 99.8) for c in (0, 1, 2)]
norm = lambda a, c: np.clip((a - lo[c]) / (hi[c] - lo[c]), 0, 1)
S = 660                                              # square crop (px), placed where the mask covers most nuclei
TS, TX, TY = 26.0, (1.0, 28.5), (8.0, 35.5)          # tile size and positions (mm)
fig.text((TX[0] + TS / 2 - 3.6) / W, 1 - 5.0 / H, 'Su(Hw)', fontsize=FS, ha='center', va='top', fontweight='bold', color='#c51b7d')
fig.text((TX[0] + TS / 2 + 3.9) / W, 1 - 5.0 / H, 'γH2Av', fontsize=FS, ha='center', va='top', fontweight='bold', color='#1a9641')
fig.text((TX[1] + TS / 2) / W, 1 - 5.0 / H, 'DAPI + nuclear mask', fontsize=FS, ha='center', va='top', fontweight='bold')
for i, ((key, lab, r), raw, mk) in enumerate(zip(EX, raws, masks)):
    cov = ndimage.uniform_filter(mk.astype(float), S, mode='constant')
    h2 = S // 2; cy, cx = np.unravel_index(np.argmax(cov[h2:-h2, h2:-h2]), cov[h2:-h2, h2:-h2].shape)
    sl = (slice(cy, cy + S), slice(cx, cx + S))
    merge = np.stack([norm(raw[..., 0][sl], 0), norm(raw[..., 1][sl], 1), norm(raw[..., 0][sl], 0)], -1)
    dapi = np.repeat(norm(raw[..., 2][sl], 2)[..., None], 3, -1)
    m = mk[sl]; edge = m & ~ndimage.binary_erosion(m, iterations=2)
    dapi[edge] = [1.0, 0.83, 0.0]                    # mask outline only on the channel it was derived from
    for j, img in enumerate((merge, dapi)):
        ia = ax_mm(TX[j], TY[i], TS, TS)
        ia.imshow(img, interpolation='nearest'); ia.set_axis_off()
        if j == 0:
            ia.text(0.03, 0.97, f'{lab}   R = {r:+.2f}', transform=ia.transAxes, ha='left', va='top', color='w', fontsize=FS,
                    fontweight='bold', bbox=dict(boxstyle='square,pad=0.25', fc='k', ec='none', alpha=0.6))
        if i == 0 and j == 0:
            bar = 10 / PX_UM; ia.plot([S - 30 - bar, S - 30], [S - 30, S - 30], color='w', lw=1.2, solid_capstyle='butt')
            ia.text(S - 30 - bar / 2, S - 48, '10 µm', color='w', ha='center', va='bottom', fontsize=FS,
                    bbox=dict(boxstyle='square,pad=0.15', fc='k', ec='none', alpha=0.6))
fig.text(1 / W, 1 - 66.2 / H, 'Yellow outline: nuclear mask Agentic-J derived from DAPI', color='#444444', fontsize=FS, va='top')

# ── b: per-field values, published vs Agentic-J ──────────────────────────────
letter(60, 0, 'b'); title(65, 0.4, "Pearson's R above Costes threshold, per field")
ax = ax_mm(71, 8, 46, 52)
rng = np.random.default_rng(3)
pos = {('iso', 'p'): 0.0, ('iso', 'a'): 0.42, ('hyp', 'p'): 1.25, ('hyp', 'a'): 1.67}
for ci, cond in enumerate(('iso', 'hyp')):
    for src, vals, col in (('p', paper[ci], COL['published']), ('a', agent[ci], COL['hil'])):
        x0 = pos[(cond, src)]
        ax.scatter(x0 + rng.uniform(-0.09, 0.09, len(vals)), vals, s=8, color=col, edgecolor='k', linewidth=0.3, zorder=3)
        mval, sd = vals.mean(), vals.std(ddof=1)
        ax.plot([x0 - 0.15, x0 + 0.15], [mval, mval], color='k', lw=0.9, zorder=4)
        ax.plot([x0, x0], [mval - sd, mval + sd], color='k', lw=0.6, zorder=4)
ax.axhline(0, color='#888888', lw=0.5, ls='--', zorder=1)
ax.set_xticks([0.21, 1.46]); ax.set_xticklabels(['Isotonic\n(n = 6 fields)', '250 mM NaCl\n(n = 11 fields)'])
ax.set_xlim(-0.35, 2.0); ax.set_ylim(-0.45, 1.0); ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_ylabel("Pearson's R (above Costes threshold)")
ax.scatter([], [], s=8, color=COL['published'], edgecolor='k', linewidth=0.3, label='Published (Fig. 5)')
ax.scatter([], [], s=8, color=COL['hil'], edgecolor='k', linewidth=0.3, label='Agentic-J (Coloc 2)')
ax.legend(loc='upper left', fontsize=FS, frameon=False, handletextpad=0.2, labelspacing=0.3, borderaxespad=0.2)

# ── c: rank-matched field values + effect summary ────────────────────────────
letter(124, 0, 'c'); title(129, 0.4, 'Field values, rank-matched')
bx = ax_mm(136, 8, 34, 34)
lim = (-0.45, 0.85)
bx.plot(lim, lim, color=GREY, lw=0.6, ls='--', zorder=1)
for ci, (mk_, lab) in enumerate((('o', 'Isotonic'), ('s', '250 mM NaCl'))):
    bx.scatter(np.sort(paper[ci]), np.sort(agent[ci]), s=9, marker=mk_, color=COL['hil'], edgecolor='k', linewidth=0.3,
               zorder=3, label=lab)
bx.set_xlim(lim); bx.set_ylim(lim); bx.set_aspect('equal')
bx.set_xticks([-0.4, 0, 0.4, 0.8]); bx.set_yticks([-0.4, 0, 0.4, 0.8])
bx.set_xlabel('Published R (sorted)'); bx.set_ylabel('Agentic-J R (sorted)')
bx.legend(loc='upper left', fontsize=FS, frameon=False, handletextpad=0.1, labelspacing=0.3, borderaxespad=0.1)
mad = np.mean(np.abs(np.concatenate([np.sort(paper[i]) - np.sort(agent[i]) for i in (0, 1)])))
bx.text(0.97, 0.04, f'mean |Δ| = {mad:.3f}', transform=bx.transAxes, ha='right', va='bottom', fontsize=FS)

lines = []
for lab, (iso, hyp), col in (('Published', paper, COL['published']), ('Agentic-J', agent, COL['hil'])):
    d = hyp.mean() - iso.mean(); ci = boot_ci(iso, hyp); p = perm_p(iso, hyp)
    lines.append((lab, col, d, ci, p, iso, hyp))
y = 52.5
fig.text(129 / W, 1 - y / H, 'Increase under NaCl (Δ mean R)', fontsize=FS, fontweight='bold', va='top')
for j, (lab, col, d, ci, p, iso, hyp) in enumerate(lines):
    fig.text(129 / W, 1 - (y + 3.2 + j * 3.0) / H, f'{lab}: {d:+.3f}  [{ci[0]:+.2f}, {ci[1]:+.2f}],  P = {fmt_p(p)}',
             fontsize=FS, color=col, va='top')
fig.text(129 / W, 1 - (y + 10.2) / H, '95% bootstrap CI; exact permutation test.\nPublished values are not linked to files,\nso fields are matched by rank.',
         fontsize=FS, color='#555555', va='top', linespacing=1.1)

for ext in ('pdf', 'svg'):
    fig.savefig(os.path.join(HERE, f'Fig_Coloc2_v2.{ext}'))
fig.savefig(os.path.join(HERE, 'Fig_Coloc2_v2.png'), dpi=600)
with open(os.path.join(HERE, 'figure_numbers_v2.txt'), 'w') as f:
    for lab, col, d, ci, p, iso, hyp in lines:
        f.write(f'{lab}: iso {iso.mean():+.3f}±{iso.std(ddof=1):.3f} (median {np.median(iso):+.3f}), hyp {hyp.mean():.3f}±{hyp.std(ddof=1):.3f} '
                f'(median {np.median(hyp):.3f}), delta {d:+.3f} (95% CI {ci[0]:+.3f} to {ci[1]:+.3f}), P = {p:.2g}\n')
    f.write(f'Rank-matched mean |difference|: {mad:.3f}\n')
    f.write(f'Example fields: isotonic field 6 (R = {EX[0][2]:+.2f}), NaCl field 1 (R = {EX[1][2]:+.2f}); pixel size {PX_UM} um\n')
print(open(os.path.join(HERE, 'figure_numbers_v2.txt')).read())
