"""BoneJ reproduction figure: Agentic-J vs Smith et al. 2023 (Biol J Linn Soc 140:1-25), Table 2.

Style as Fig. 2: Nature sizes (min 5 pt), 180 mm wide, published = purple, Agentic-J (automated) = green.
Inputs (inputs/): reference_table2.csv (paper species means + CV + body mass), species_summary.csv and
specimens.csv (Agentic-J run 2026-10-05, instance 1), tbsp_maskArtefacts_test.tsv (BoneJ re-run of two
species with maskArtefacts on/off). Outputs: Fig_BoneJ_reproduction.{pdf,svg,png} next to this script.
"""
import csv, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); IN = os.path.join(HERE, 'inputs')
FS, FL = 5, 6
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': FS,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'xtick.labelsize': FS, 'ytick.labelsize': FS, 'axes.labelsize': FS, 'axes.titlesize': FS, 'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5})
COL = {'published': '#7b3294', 'automated': '#1a9641'}
LIGHT = '#a6dba0'                     # per-specimen Agentic-J values
GREY = '#bbbbbb'
MM = 1 / 25.4; W, H = 180, 150
fig = plt.figure(figsize=(W * MM, H * MM))
def ax_mm(x, y, w, h, **k): return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **k)
def letter(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight='bold', va='top')
def title(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=FL, va='top')

# ── data ─────────────────────────────────────────────────────────────────────
ref = list(csv.DictReader(open(f'{IN}/reference_table2.csv')))
ag = {r['species'].split()[-1]: r for r in csv.DictReader(open(f'{IN}/species_summary.csv'))}
spec = list(csv.DictReader(open(f'{IN}/specimens.csv')))
mask = list(csv.DictReader(open(f'{IN}/tbsp_maskArtefacts_test.tsv'), delimiter='\t'))
METRICS = [  # key in Table 2, column stem in Agentic-J tables, label, unit, axis range, ticks
    ('bv_tv', 'bone_volume_fraction', 'BV/TV', '', (0.16, 0.36), [0.2, 0.25, 0.3, 0.35]),
    ('tb_th_mm', 'mean_trabecular_thickness_mm', 'Tb.Th', 'mm', (0.06, 0.14), [0.06, 0.08, 0.10, 0.12, 0.14]),
    ('tb_sp_mm', 'mean_trabecular_spacing_mm', 'Tb.Sp', 'mm', (0.22, 0.44), [0.25, 0.30, 0.35, 0.40]),
    ('conn_d_per_mm3', 'connectivity_density_per_mm3', 'Conn.D', 'mm$^{-3}$', (0, 80), [0, 20, 40, 60, 80]),
]
CV = {'bv_tv': 'bv_tv_cv', 'tb_th_mm': 'tb_th_cv', 'tb_sp_mm': 'tb_sp_cv', 'conn_d_per_mm3': 'conn_d_cv'}
species = [r['species'] for r in ref]
mass = np.array([float(r['mass_g']) for r in ref])
P = {m: np.array([float(r[m]) for r in ref]) for m, *_ in METRICS}
Psd = {m: P[m] * np.array([float(r[CV[m]]) for r in ref]) for m, *_ in METRICS}
A = {m: np.array([float(ag[s][k + '_mean']) for s in species]) for m, k, *_ in METRICS}
Asd = {m: np.array([float(ag[s][k + '_sd']) for s in species]) for m, k, *_ in METRICS}
massof = {r['genus'] + ' ' + r['species']: float(r['mass_g']) for r in ref}

def lin_ccc(x, y):
    mx, my = x.mean(), y.mean(); sxy = ((x - mx) * (y - my)).mean()
    return 2 * sxy / (x.var() + y.var() + (mx - my) ** 2)

# masking-off re-run: species means per species
mask_off = {}
for sp in {r['species'] for r in mask}:
    v = [float(r['mask_false']) for r in mask if r['species'] == sp]
    mask_off[sp.split()[-1]] = np.mean(v)

def style(ax, unit_label, xlab=None):
    ax.tick_params(length=2, width=0.5)
    if xlab: ax.set_xlabel(xlab)

# ── a: agreement of species means ────────────────────────────────────────────
letter(0, 0, 'a'); title(5, 0.4, 'Species means, Agentic-J vs published (11 species, Table 2)')
pw, gap, x0, y0, ph = 35.0, 10.0, 8.0, 7.0, 34.0
for i, (m, k, lab, unit, lim, ticks) in enumerate(METRICS):
    ax = ax_mm(x0 + i * (pw + gap), y0, pw, ph)
    ax.plot(lim, lim, color=GREY, lw=0.6, ls='--', zorder=1)
    ax.errorbar(P[m], A[m], xerr=Psd[m], yerr=Asd[m], fmt='none', ecolor='#d9d9d9', elinewidth=0.4, capsize=0, zorder=2)
    ax.scatter(P[m], A[m], s=9, color=COL['automated'], edgecolor='k', linewidth=0.3, zorder=3, label='As run')
    if m == 'tb_sp_mm':
        for sp, v in mask_off.items():
            j = species.index(sp)
            ax.scatter(P[m][j], v, s=11, facecolor='w', edgecolor=COL['automated'], linewidth=0.7, zorder=4, marker='D',
                       label='Masking off (re-run)' if sp == sorted(mask_off)[0] else None)
            ax.annotate('', xy=(P[m][j], v), xytext=(P[m][j], A[m][j]), zorder=3,
                        arrowprops=dict(arrowstyle='-|>', color=COL['automated'], lw=0.5, mutation_scale=4, shrinkA=2, shrinkB=2.5))
        ax.legend(loc='upper left', fontsize=FS, frameon=False, handletextpad=0.2, borderaxespad=0.2, labelspacing=0.3)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xticks(ticks); ax.set_yticks(ticks); ax.set_aspect('equal')
    u = f' ({unit})' if unit else ''
    ax.set_xlabel(f'Published {lab}{u}'); ax.set_ylabel(f'Agentic-J {lab}{u}')
    d = 100 * (A[m] - P[m]) / P[m]
    txt = f'CCC = {lin_ccc(P[m], A[m]):.4f}\nmean |Δ| = {np.mean(np.abs(d)):.1f}%'
    ax.text(0.97, 0.04, txt, transform=ax.transAxes, ha='right', va='bottom', fontsize=FS, linespacing=1.15)
    ax.set_title(lab, fontsize=FL, pad=2)

# ── b: relative deviation per species ────────────────────────────────────────
YB = 52
letter(0, YB, 'b'); title(5, YB + 0.4, 'Deviation from published species means')
ax = ax_mm(16, YB + 8, 74, 32)
ax.axvspan(-1, 1, color='#eeeeee', lw=0, zorder=0)
ax.axvline(0, color='#888888', lw=0.5, zorder=1)
rng = np.random.default_rng(1)
for r, (m, k, lab, *_) in enumerate(METRICS):
    d = 100 * (A[m] - P[m]) / P[m]
    jit = rng.uniform(-0.17, 0.17, len(d))
    ax.scatter(d, r + jit, s=8, color=COL['automated'], edgecolor='k', linewidth=0.3, zorder=3)
    if m == 'tb_sp_mm':
        for sp, v in mask_off.items():
            j = species.index(sp)
            ax.scatter(100 * (v - P[m][j]) / P[m][j], r - 0.34, s=10, facecolor='w', edgecolor=COL['automated'],
                       linewidth=0.7, marker='D', zorder=4)
        ax.text(1.4, r - 0.34, 'masking off (re-run, 2 species)', va='center', ha='left', fontsize=FS, color='#444444')
ax.set_yticks(range(len(METRICS))); ax.set_yticklabels([m[2] for m in METRICS]); ax.invert_yaxis()
ax.set_xlim(-2, 12); ax.set_xticks([-2, 0, 2, 4, 6, 8, 10, 12])
ax.set_xlabel('Agentic-J − published (% of published mean)')
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
ax.set_ylim(len(METRICS) - 0.5, -0.75)
ax.text(1.1, -0.62, '±1% band', ha='left', va='center', fontsize=FS, color='#777777')

# ── c: allometric scaling exponents (the biological result) ──────────────────
letter(98, YB, 'c'); title(103, YB + 0.4, 'Scaling with body mass, log–log slope ± 95% CI')
ax = ax_mm(112, YB + 8, 66, 32)
lm = np.log10(mass)
def slope_ci(y):
    res = stats.linregress(lm, np.log10(y)); t = stats.t.ppf(0.975, len(y) - 2)
    return res.slope, t * res.stderr, res.pvalue
for r, (m, k, lab, *_) in enumerate(METRICS):
    for off, (name, src) in zip((-0.14, 0.14), (('published', P[m]), ('automated', A[m]))):
        s, ci, p = slope_ci(src)
        ax.errorbar(s, r + off, xerr=ci, fmt='o', ms=2.6, color=COL[name], mec='k', mew=0.3, elinewidth=0.7, capsize=1.2, capthick=0.5)
ax.axvline(0, color='#888888', lw=0.5)
ax.set_yticks(range(len(METRICS))); ax.set_yticklabels([m[2] for m in METRICS]); ax.tick_params(axis='y', length=0)
ax.spines['left'].set_visible(False)
ax.set_xlim(-1.1, 1.0); ax.set_xticks([-1, -0.5, 0, 0.5, 1.0])
ax.set_xlabel('Slope of log$_{10}$(metric) on log$_{10}$(body mass)')
ax.scatter([], [], s=8, color=COL['published'], edgecolor='k', linewidth=0.3, label='Published (Table 2)')
ax.scatter([], [], s=8, color=COL['automated'], edgecolor='k', linewidth=0.3, label='Agentic-J')
ax.legend(loc='upper right', fontsize=FS, frameon=False, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.3)
ax.set_ylim(len(METRICS) - 0.5, -0.6)

# ── d: metric vs body mass, per specimen + species means ─────────────────────
YD = 99
letter(0, YD, 'd'); title(5, YD + 0.4, 'Trabecular architecture vs body mass (66 specimens, 11 species)')
xs = np.array([20, 200])
for i, (m, k, lab, unit, lim, ticks) in enumerate(METRICS):
    ax = ax_mm(x0 + i * (pw + gap), YD + 9, pw, 34)
    xm = np.array([massof[r['species']] for r in spec]); yv = np.array([float(r[k]) for r in spec])
    jit = 10 ** rng.uniform(-0.012, 0.012, len(xm))
    ax.scatter(xm * jit, yv, s=3, color=LIGHT, lw=0, zorder=1, label='Agentic-J, per specimen')
    for name, src in (('published', P[m]), ('automated', A[m])):
        res = stats.linregress(np.log10(mass), np.log10(src))
        ax.plot(xs, 10 ** (res.intercept + res.slope * np.log10(xs)), color=COL[name], lw=0.8,
                ls='-' if name == 'published' else (0, (3, 1.5)), zorder=2)
    ax.scatter(mass, P[m], s=12, facecolor='none', edgecolor=COL['published'], linewidth=0.7, zorder=3, label='Published, species mean')
    ax.scatter(mass, A[m], s=5, color=COL['automated'], edgecolor='k', linewidth=0.25, zorder=4, label='Agentic-J, species mean')
    ax.set_xscale('log'); ax.set_xlim(20, 200); ax.set_xticks([20, 50, 100, 200]); ax.set_xticklabels(['20', '50', '100', '200'])
    ax.minorticks_off()
    ax.set_xlabel('Body mass (g)'); ax.set_ylabel(f'{lab}' + (f' ({unit})' if unit else ''))
    if i == 0:
        h_, l_ = ax.get_legend_handles_labels()
        from matplotlib.lines import Line2D
        h_ += [Line2D([], [], color=COL['published'], lw=0.8), Line2D([], [], color=COL['automated'], lw=0.8, ls=(0, (3, 1.5)))]
        l_ += ['Published fit', 'Agentic-J fit']
        fig.legend(h_, l_, loc='upper right', bbox_to_anchor=(1, 1 - (YD + 0.2) / H), ncol=5, fontsize=FS, frameon=False,
                   handletextpad=0.3, columnspacing=0.9, handlelength=1.6, borderaxespad=0.2)

for ext in ('pdf', 'svg'):
    fig.savefig(os.path.join(HERE, f'Fig_BoneJ_reproduction.{ext}'))
fig.savefig(os.path.join(HERE, 'Fig_BoneJ_reproduction.png'), dpi=600)

# numbers quoted in the legend
with open(os.path.join(HERE, 'figure_numbers.txt'), 'w') as f:
    for m, k, lab, *_ in METRICS:
        d = 100 * (A[m] - P[m]) / P[m]
        sp_, cp_, pp_ = slope_ci(P[m]); sa_, ca_, pa_ = slope_ci(A[m])
        f.write(f'{lab}: CCC {lin_ccc(P[m], A[m]):.4f}; mean|d| {np.mean(np.abs(d)):.2f}%; max|d| {np.max(np.abs(d)):.2f}%; '
                f'slope published {sp_:.3f}±{cp_:.3f} (p={pp_:.2g}), Agentic-J {sa_:.3f}±{ca_:.3f} (p={pa_:.2g})\n')
    for sp, v in mask_off.items():
        j = species.index(sp)
        f.write(f'Tb.Sp masking off {sp}: {v:.4f} vs published {P["tb_sp_mm"][j]:.3f} ({100*(v-P["tb_sp_mm"][j])/P["tb_sp_mm"][j]:+.2f}%)\n')
print(open(os.path.join(HERE, 'figure_numbers.txt')).read())
