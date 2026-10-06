"""Coloc 2 reproduction figure: Agentic-J vs Simmons et al. 2022 (PLoS Genet, doi:10.1371/journal.pgen.1010396), Fig. 5A-C.

Style as Fig. 2: Nature sizes (min 5 pt), 180 mm wide; published = purple, interactive (human-in-the-loop) = blue,
auto-pilot = green, fast mode = brown. Inputs (inputs/):
  paper_fig5_pcc.csv            Pearson's R above Costes threshold per field, from the authors' 'Figure 5 data.xlsx'
  v2_interactive_coloc2.csv     Agentic-J, prompt v2, interactive, Coloc 2 (2026-10-06)
  v1_interactive_coloc2.csv     Agentic-J, prompt v1, interactive, Coloc 2 (2026-10-05; above-threshold R read from its Coloc 2 log)
  v1_advanced_autopilot.csv     Agentic-J, prompt v1, auto-pilot, cp_measure (no-threshold Pearson only)
  v1_quick_autopilot.csv        Agentic-J, prompt v1, fast mode, custom Groovy (no-threshold Pearson only)
Outputs: Fig_Coloc2_reproduction.{pdf,svg,png} and figure_numbers.txt next to this script.
"""
import csv, itertools, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); IN = os.path.join(HERE, 'inputs')
FS, FL = 5, 6
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': FS,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'xtick.labelsize': FS, 'ytick.labelsize': FS, 'axes.labelsize': FS, 'axes.titlesize': FS, 'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'xtick.major.pad': 1.2, 'ytick.major.pad': 1.2, 'axes.labelpad': 1.5})
COL = {'published': '#7b3294', 'hil': '#0571b0', 'automated': '#1a9641', 'fast': '#a6611a'}
GREY = '#bbbbbb'
MM = 1 / 25.4; W, H = 180, 72
fig = plt.figure(figsize=(W * MM, H * MM))
def ax_mm(x, y, w, h, **k): return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **k)
def letter(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=8, fontweight='bold', va='top')
def title(x, y, s): fig.text(x / W, 1 - y / H, s, fontsize=FL, va='top')

# ── data ─────────────────────────────────────────────────────────────────────
def load(name, col, cond_col='condition'):
    rows = list(csv.DictReader(open(f'{IN}/{name}')))
    iso = np.array([float(r[col]) for r in rows if r[cond_col].startswith('iso')])
    hyp = np.array([float(r[col]) for r in rows if r[cond_col].startswith('hyper')])
    return iso, hyp
paper = load('paper_fig5_pcc.csv', 'pearson_r_above_threshold')
v2_above = load('v2_interactive_coloc2.csv', 'pearson_r_above_costes')
v2_all = load('v2_interactive_coloc2.csv', 'pearson_r_all_nuclear')
v1_above = load('v1_interactive_coloc2.csv', 'pearson_above_threshold')
v1_all = load('v1_interactive_coloc2.csv', 'pearson_no_threshold')
auto_all = load('v1_advanced_autopilot.csv', 'pearson_r')
fast_all = load('v1_quick_autopilot.csv', 'pearson_r')

def perm_p(x, y):
    v = np.concatenate([x, y]); obs = abs(y.mean() - x.mean()); n = tot = 0
    for c in itertools.combinations(range(len(v)), len(x)):
        m = np.zeros(len(v), bool); m[list(c)] = True
        n += abs(v[~m].mean() - v[m].mean()) >= obs - 1e-12; tot += 1
    return n / tot
def fmt_p(p):
    m, e = f'{p:.1e}'.split('e'); return f'{m}×10$^{{{int(e)}}}$'
def boot_ci(x, y, n=10000, seed=20261006):
    rng = np.random.default_rng(seed)
    d = rng.choice(y, (n, len(y))).mean(1) - rng.choice(x, (n, len(x))).mean(1)
    return np.percentile(d, [2.5, 97.5])

# ── a: per-field values, published vs Agentic-J (prompt v2) ──────────────────
letter(0, 0, 'a'); title(5, 0.4, "Pearson's R above Costes threshold, per image field")
ax = ax_mm(11, 9, 52, 52)
rng = np.random.default_rng(3)
pos = {('iso', 'p'): 0.0, ('iso', 'a'): 0.42, ('hyp', 'p'): 1.25, ('hyp', 'a'): 1.67}
for ci, cond in enumerate(('iso', 'hyp')):
    for src, vals, col in (('p', paper[ci], COL['published']), ('a', v2_above[ci], COL['hil'])):
        x0 = pos[(cond, src)]
        ax.scatter(x0 + rng.uniform(-0.09, 0.09, len(vals)), vals, s=8, color=col, edgecolor='k', linewidth=0.3, zorder=3)
        m, sd = vals.mean(), vals.std(ddof=1)
        ax.plot([x0 - 0.15, x0 + 0.15], [m, m], color='k', lw=0.9, zorder=4)
        ax.plot([x0, x0], [m - sd, m + sd], color='k', lw=0.6, zorder=4)
ax.axhline(0, color='#888888', lw=0.5, ls='--', zorder=1)
ax.set_xticks([0.21, 1.46]); ax.set_xticklabels(['Isotonic\n(n = 6 fields)', '250 mM NaCl\n(n = 11 fields)'])
ax.set_xlim(-0.35, 2.0); ax.set_ylim(-0.45, 1.0); ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_ylabel("Pearson's R (above Costes threshold)")
for src, (iso, hyp), col, y in (('Published', paper, COL['published'], 0.95), ('Agentic-J', v2_above, COL['hil'], 0.86)):
    ax.text(1.98, y, f'{src}: Δ = {hyp.mean() - iso.mean():+.3f}, P = {fmt_p(perm_p(iso, hyp))}', ha='right', va='center',
            fontsize=FS, color=col)
ax.scatter([], [], s=8, color=COL['published'], edgecolor='k', linewidth=0.3, label='Published (Fig. 5)')
ax.scatter([], [], s=8, color=COL['hil'], edgecolor='k', linewidth=0.3, label='Agentic-J (Coloc 2)')
ax.legend(loc='upper left', bbox_to_anchor=(0, 0.83), fontsize=FS, frameon=False, handletextpad=0.2, labelspacing=0.3, borderaxespad=0.2)

# ── b: quantile-quantile, Agentic-J vs published ─────────────────────────────
letter(68, 0, 'b'); title(73, 0.4, 'Field values, rank-matched')
bx = ax_mm(78, 9, 34, 34)
lim = (-0.45, 0.85)
bx.plot(lim, lim, color=GREY, lw=0.6, ls='--', zorder=1)
for ci, (mk, lab) in enumerate((('o', 'Isotonic'), ('s', '250 mM NaCl'))):
    p, a = np.sort(paper[ci]), np.sort(v2_above[ci])
    bx.scatter(p, a, s=9, marker=mk, color=COL['hil'], edgecolor='k', linewidth=0.3, zorder=3, label=lab)
bx.set_xlim(lim); bx.set_ylim(lim); bx.set_aspect('equal')
bx.set_xticks([-0.4, 0, 0.4, 0.8]); bx.set_yticks([-0.4, 0, 0.4, 0.8])
bx.set_xlabel("Published R (sorted)"); bx.set_ylabel("Agentic-J R (sorted)")
bx.legend(loc='upper left', fontsize=FS, frameon=False, handletextpad=0.1, labelspacing=0.3, borderaxespad=0.1)
mad = np.mean(np.abs(np.concatenate([np.sort(paper[i]) - np.sort(v2_above[i]) for i in (0, 1)])))
bx.text(0.97, 0.04, f'mean |Δ| = {mad:.3f}', transform=bx.transAxes, ha='right', va='bottom', fontsize=FS)
fig.text(78 / W, 1 - 52.5 / H, 'The published sheet does not link values\nto files; fields are matched by rank.',
         fontsize=FS, va='top', color='#555555', linespacing=1.1)

# ── c: effect of osmotic stress across runs and metrics ──────────────────────
letter(118, 0, 'c'); title(123, 0.4, 'Effect of osmotic stress (Δ mean R, 95% CI)')
cx = ax_mm(150, 9, 28, 52)
RUNS = [  # label, data, colour, filled (= published metric)
    ('Published (Coloc 2)', paper, COL['published'], True),
    ('Interactive, prompt v2 (Coloc 2)', v2_above, COL['hil'], True),
    ('Interactive, prompt v1 (Coloc 2)', v1_above, COL['hil'], True),
    (None, None, None, None),
    ('Interactive, prompt v2 (Coloc 2)', v2_all, COL['hil'], False),
    ('Interactive, prompt v1 (Coloc 2)', v1_all, COL['hil'], False),
    ('Auto-pilot (cp_measure)', auto_all, COL['automated'], False),
    ('Fast mode (custom Groovy)', fast_all, COL['fast'], False),
]
numbers = []
for r, (lab, data, col, filled) in enumerate(RUNS):
    if lab is None:
        cx.axhline(r, color='#dddddd', lw=0.5); continue
    iso, hyp = data; d = hyp.mean() - iso.mean(); lo, hi = boot_ci(iso, hyp); p = perm_p(iso, hyp)
    cx.errorbar(d, r, xerr=[[d - lo], [hi - d]], fmt='o', ms=3, mfc=col if filled else 'w', mec=col, mew=0.8,
                ecolor=col, elinewidth=0.7, capsize=1.2, capthick=0.5)
    cx.text(-0.03, r, lab, transform=cx.get_yaxis_transform(), ha='right', va='center', fontsize=FS)
    numbers.append(f'{lab} [{"above threshold" if filled else "no threshold"}]: iso {iso.mean():+.3f}±{iso.std(ddof=1):.3f}, '
                   f'hyp {hyp.mean():.3f}±{hyp.std(ddof=1):.3f}, Δ {d:+.3f} (95% CI {lo:+.3f} to {hi:+.3f}), P = {p:.2g}')
cx.text(-0.03, -0.75, 'Above Costes threshold (published metric)', transform=cx.get_yaxis_transform(), ha='right', va='center',
        fontsize=FS, fontweight='bold')
cx.text(-0.03, 3.4, 'Without threshold', transform=cx.get_yaxis_transform(), ha='right', va='center', fontsize=FS, fontweight='bold')
pub_d = paper[1].mean() - paper[0].mean()
cx.axvline(pub_d, color=COL['published'], lw=0.5, ls=':', zorder=0)
cx.axvline(0, color='#888888', lw=0.5)
cx.set_ylim(len(RUNS) - 0.5, -1.2); cx.set_yticks([])
cx.spines['left'].set_visible(False)
cx.set_xlim(-0.05, 0.8); cx.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
cx.set_xlabel('NaCl − isotonic')

for ext in ('pdf', 'svg'):
    fig.savefig(os.path.join(HERE, f'Fig_Coloc2_reproduction.{ext}'))
fig.savefig(os.path.join(HERE, 'Fig_Coloc2_reproduction.png'), dpi=600)
numbers.append(f'Panel b rank-matched mean |difference| (prompt v2 vs published): {mad:.3f}')
open(os.path.join(HERE, 'figure_numbers.txt'), 'w').write('\n'.join(numbers) + '\n')
print('\n'.join(numbers))
