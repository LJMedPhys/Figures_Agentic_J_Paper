# Stage 2 (Plotting only)
# Publication-quality single-category beeswarm plot of per-interval mean velocities.
# Inputs: interval_mean_velocities.csv (per-interval mean velocities)
# Outputs: figures/interval_mean_velocity_beeswarm_padded_ylabel28_yticks_big_dots_ylim3_size530x600.(png|svg)

# NOTE: Environment presets already provide: pd, np, plt, sns, stats, os

# ---------------------------
# User parameters / paths
# ---------------------------
PROJECT_DIR = r"/app/data/projects/mouse_cell_tracking_trackmate"
INPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv"
FIG_DIR = os.path.join(PROJECT_DIR, "figures")

OUT_BASENAME = "interval_mean_velocity_beeswarm_padded_ylabel28_yticks_big_dots_ylim3_size530x600"
OUT_PNG = os.path.join(FIG_DIR, f"{OUT_BASENAME}.png")
OUT_SVG = os.path.join(FIG_DIR, f"{OUT_BASENAME}.svg")

# Exact output pixel size requirement
TARGET_W_PX = 530
TARGET_H_PX = 600
DPI = 300
FIG_W_IN = TARGET_W_PX / DPI
FIG_H_IN = TARGET_H_PX / DPI

# Styling parameters
Y_LABEL = "Mean velocity (µm/min)"
Y_LABEL_FONTSIZE = 28
Y_TICK_FONTSIZE = 22

# Larger dots while keeping swarm packing
DOT_SIZE = 10  # seaborn swarmplot 'size' is in points
DOT_COLOR = "#0173B2"  # colorblind-safe blue

# Mean segment styling
MEAN_LINEWIDTH = 3.0
MEAN_COLOR = "black"

# Axis limits/ticks
Y_LIM = (0, 3)
Y_TICKS = [0, 1, 2, 3]

# Padding / layout
LEFT_PAD = 0.22
RIGHT_PAD = 0.98
BOTTOM_PAD = 0.08
TOP_PAD = 0.98

# ---------------------------
# Load data
# ---------------------------
os.makedirs(FIG_DIR, exist_ok=True)

df = pd.read_csv(INPUT_CSV)

required_cols = {"mean_velocity_um_per_min"}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in {INPUT_CSV}: {missing}. Found columns: {list(df.columns)}")

vals = pd.to_numeric(df["mean_velocity_um_per_min"], errors="coerce").dropna()
if vals.empty:
    raise ValueError("No valid numeric values found in mean_velocity_um_per_min after coercion.")

overall_mean = float(vals.mean())

print(f"N intervals (rows) used for plot: {len(vals)}")
print(f"Overall mean of interval mean velocities (µm/min): {overall_mean:.6f}")

# ---------------------------
# Plot
# ---------------------------
sns.set_style("white")
sns.set_palette("colorblind")

plt.rcParams.update({
    "font.size": 12,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
})

fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN), dpi=DPI)

plot_df = pd.DataFrame({
    "group": [""] * len(vals),
    "mean_velocity_um_per_min": vals.values
})

# Single-category beeswarm
sns.swarmplot(
    data=plot_df,
    x="group",
    y="mean_velocity_um_per_min",
    ax=ax,
    color=DOT_COLOR,
    size=DOT_SIZE,
    edgecolor=None,
    linewidth=0,
)

# Short overall-mean segment spanning the dot column
# For a single category at x=0, draw a short segment centered at 0.
segment_half_width = 0.12
ax.plot(
    [-segment_half_width, segment_half_width],
    [overall_mean, overall_mean],
    color=MEAN_COLOR,
    linewidth=MEAN_LINEWIDTH,
    solid_capstyle="butt",
    zorder=5,
)

# Axes formatting
ax.set_xlabel("")
ax.set_ylabel(Y_LABEL, fontsize=Y_LABEL_FONTSIZE)

# No x ticks/labels
ax.set_xticks([])
ax.set_xticklabels([])

# Y axis limits and ticks
ax.set_ylim(*Y_LIM)
ax.set_yticks(Y_TICKS)

# Larger y tick labels + visible tick marks on left, pointing outward
ax.tick_params(
    axis="y",
    which="major",
    labelsize=Y_TICK_FONTSIZE,
    left=True,
    right=False,
    direction="out",
    length=6,
    width=1.5,
)

# Remove grid
ax.grid(False)

# Minimal spines: left/bottom only
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# Add left/right padding to avoid crowding axes
# For a single category at x=0, set xlim with padding.
ax.set_xlim(-0.5, 0.5)

# Tight layout with explicit padding
fig.subplots_adjust(left=LEFT_PAD, right=RIGHT_PAD, bottom=BOTTOM_PAD, top=TOP_PAD)

# Save with exact pixel dimensions (figsize*dpi)
fig.savefig(OUT_PNG, dpi=DPI, bbox_inches=None, pad_inches=0)
fig.savefig(OUT_SVG, bbox_inches=None, pad_inches=0)
plt.close(fig)

# Verify pixel dimensions for PNG
try:
    from PIL import Image
    with Image.open(OUT_PNG) as im:
        w, h = im.size
    print(f"Saved PNG pixel size: {w}x{h} (target {TARGET_W_PX}x{TARGET_H_PX})")
except Exception as e:
    print(f"PNG size verification skipped/failed (PIL not available?): {e}")

print(f"Saved: {OUT_PNG}")
print(f"Saved: {OUT_SVG}")
