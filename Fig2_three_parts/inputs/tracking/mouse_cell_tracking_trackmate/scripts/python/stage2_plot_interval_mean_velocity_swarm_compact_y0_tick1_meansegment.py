# Stage 2 (plotting only): Compact single-category swarm/beeswarm plot of per-interval mean velocities
# Input:  /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv
# Output: /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_swarm_compact_y0_tick1_meansegment.(png|svg)

# NOTE: Do NOT run statistical tests here; this script only visualizes values.

# Environment presets (do not re-import): pd, np, plt, sns, stats, os

# -----------------------------
# Parameters
# -----------------------------
INPUT_CSV = r'/app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv'
PROJECT_DIR = r'/app/data/projects/mouse_cell_tracking_trackmate'
FIG_DIR = os.path.join(PROJECT_DIR, 'figures')

Y_COL = 'mean_velocity_um_per_min'

# Styling parameters
POINT_COLOR = '#0173B2'  # colorblind-safe blue
MEAN_LINE_COLOR = 'black'
FIGSIZE = (3.0, 4.0)     # compact
DPI = 300

# Mean segment span around swarm center (single category at x=0)
MEAN_SEG_XMIN = -0.08
MEAN_SEG_XMAX = 0.08

# Y-axis ticks
Y_MIN = 0.0
Y_TICK_STEP = 1.0

# -----------------------------
# Load data
# -----------------------------
os.makedirs(FIG_DIR, exist_ok=True)

df = pd.read_csv(INPUT_CSV)

if Y_COL not in df.columns:
    raise ValueError(f"Expected column '{Y_COL}' not found. Columns: {list(df.columns)}")

# Coerce to numeric and drop missing
vals = pd.to_numeric(df[Y_COL], errors='coerce').dropna()

if len(vals) == 0:
    raise ValueError(f"No valid numeric values found in column '{Y_COL}'.")

# Report IQR outliers (do not remove; just report)
q1 = np.percentile(vals, 25)
q3 = np.percentile(vals, 75)
iqr = q3 - q1
lower = q1 - 1.5 * iqr
upper = q3 + 1.5 * iqr
n_outliers = int(((vals < lower) | (vals > upper)).sum())
print(f"IQR outlier check for {Y_COL}: Q1={q1:.4f}, Q3={q3:.4f}, IQR={iqr:.4f}, lower={lower:.4f}, upper={upper:.4f}, n_outliers={n_outliers} (not removed)")

overall_mean = float(vals.mean())
print(f"Overall mean {Y_COL} (used for mean segment): {overall_mean:.6f} µm/min")

# -----------------------------
# Plot
# -----------------------------
plt.close('all')
sns.set_style('white')
sns.set_palette('colorblind')

fig, ax = plt.subplots(figsize=FIGSIZE)

plot_df = pd.DataFrame({'category': [''] * len(vals), Y_COL: vals.values})

# Single-category beeswarm
sns.swarmplot(
    data=plot_df,
    x='category',
    y=Y_COL,
    ax=ax,
    color=POINT_COLOR,
    size=5,
    edgecolor=None,
    linewidth=0
)

# Short overall-mean segment centered on swarm
ax.plot([MEAN_SEG_XMIN, MEAN_SEG_XMAX], [overall_mean, overall_mean], color=MEAN_LINE_COLOR, linewidth=2.0, solid_capstyle='butt')

# Labels and axes formatting
ax.set_xlabel('')
ax.set_ylabel('Mean velocity per interval (µm/min)', fontsize=14)

# Remove x-axis ticks/labels entirely
ax.set_xticks([])
ax.set_xticklabels([])
ax.tick_params(axis='x', bottom=False)

# Y-axis starts at 0 with 1.0 tick spacing
ax.set_ylim(bottom=Y_MIN)

# Set y ticks at integer spacing up to max (rounded up)
ymax = float(np.nanmax(vals))
ymax_ceil = np.ceil(ymax / Y_TICK_STEP) * Y_TICK_STEP
# Ensure at least one tick above 0
if ymax_ceil < Y_TICK_STEP:
    ymax_ceil = Y_TICK_STEP
ax.set_yticks(np.arange(Y_MIN, ymax_ceil + 0.5 * Y_TICK_STEP, Y_TICK_STEP))
ax.tick_params(axis='y', labelsize=12)

# Minimal styling: no grid; only left/bottom spines
ax.grid(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.0)
ax.spines['bottom'].set_linewidth(1.0)

# Tight layout and save
png_path = os.path.join(FIG_DIR, 'interval_mean_velocity_swarm_compact_y0_tick1_meansegment.png')
svg_path = os.path.join(FIG_DIR, 'interval_mean_velocity_swarm_compact_y0_tick1_meansegment.svg')

plt.savefig(png_path, dpi=DPI, bbox_inches='tight')
plt.savefig(svg_path, bbox_inches='tight')

print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")
