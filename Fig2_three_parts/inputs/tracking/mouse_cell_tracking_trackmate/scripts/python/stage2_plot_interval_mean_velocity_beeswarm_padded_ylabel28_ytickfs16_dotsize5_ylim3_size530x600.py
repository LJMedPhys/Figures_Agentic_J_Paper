# Stage 2 (Plotting only)
# Regenerate publication-quality single-category beeswarm plot of per-interval mean velocities.
# Inputs:
#   /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv
# Outputs:
#   /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_beeswarm_padded_ylabel28_ytickfs16_dotsize5_ylim3_size530x600.png
#   /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_beeswarm_padded_ylabel28_ytickfs16_dotsize5_ylim3_size530x600.svg
# Notes:
# - No statistical testing performed here (Stage 2 only). Overall mean is computed for plotting only.

import os

# -----------------------------
# Config
# -----------------------------
PROJECT_DIR = r"/app/data/projects/mouse_cell_tracking_trackmate"
INPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv"
FIG_DIR = os.path.join(PROJECT_DIR, "figures")

OUT_PNG = os.path.join(
    FIG_DIR,
    "interval_mean_velocity_beeswarm_padded_ylabel28_ytickfs16_dotsize5_ylim3_size530x600.png",
)
OUT_SVG = os.path.join(
    FIG_DIR,
    "interval_mean_velocity_beeswarm_padded_ylabel28_ytickfs16_dotsize5_ylim3_size530x600.svg",
)

Y_COL = "mean_velocity_um_per_min"
Y_LABEL = "Mean velocity (µm/min)"

# Exact output size requirement: 530x600 px
DPI = 300
FIG_W_IN = 530 / DPI
FIG_H_IN = 600 / DPI

# Styling requirements
DOT_COLOR = "#0173B2"  # colorblind-safe blue
DOT_SIZE = 5  # EXACT marker size requested
MEAN_LINE_COLOR = "black"
MEAN_LINE_WIDTH = 2.0

YLIM = (0, 3)
YTICKS = [0, 1, 2, 3]
YLABEL_FONTSIZE = 28
YTICKLABEL_FONTSIZE = 16

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(INPUT_CSV)

required = {Y_COL}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in {INPUT_CSV}: {sorted(missing)}. Found: {list(df.columns)}")

# Coerce numeric and drop NaNs
vals = pd.to_numeric(df[Y_COL], errors="coerce")
vals = vals.dropna()
if vals.empty:
    raise ValueError(f"No valid numeric values found in column '{Y_COL}' from {INPUT_CSV}.")

overall_mean = float(vals.mean())
print(f"Overall mean of interval mean velocities ({Y_COL}) = {overall_mean:.6g} µm/min (N intervals={len(vals)})")

# -----------------------------
# Plot
# -----------------------------
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_palette("colorblind")

# Use a single category to create a beeswarm column
plot_df = pd.DataFrame({"category": [""] * len(vals), Y_COL: vals.values})

plt.close("all")
fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN), dpi=DPI)

# Swarmplot (beeswarm)
# - size is in points^2 for scatter-like markers in seaborn; seaborn's `size` maps to marker size.
#   We set it exactly as requested.
sns.swarmplot(
    data=plot_df,
    x="category",
    y=Y_COL,
    ax=ax,
    color=DOT_COLOR,
    size=DOT_SIZE,
    edgecolor=None,
    linewidth=0,
)

# Overall mean as a short segment spanning the dot column
# x positions for single category are centered at 0
x_center = 0
segment_half_width = 0.12
ax.plot(
    [x_center - segment_half_width, x_center + segment_half_width],
    [overall_mean, overall_mean],
    color=MEAN_LINE_COLOR,
    linewidth=MEAN_LINE_WIDTH,
    solid_capstyle="butt",
    zorder=5,
)

# Axes formatting
ax.set_xlabel("")
ax.set_ylabel(Y_LABEL, fontsize=YLABEL_FONTSIZE)

ax.set_ylim(*YLIM)
ax.set_yticks(YTICKS)
ax.tick_params(
    axis="y",
    labelsize=YTICKLABEL_FONTSIZE,
    direction="out",
    length=6,
    width=1.2,
    left=True,
    right=False,
)

# Remove x-axis ticks/labels entirely
ax.set_xticks([])
ax.tick_params(axis="x", bottom=False, labelbottom=False)

# Minimal spines: left + bottom only
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["left", "bottom"]:
    ax.spines[spine].set_linewidth(1.2)

# No grid
ax.grid(False)

# Padding to avoid crowding axes
# Increase left/right margins; keep tight vertical packing
plt.subplots_adjust(left=0.28, right=0.98, bottom=0.06, top=0.98)

# Save
plt.savefig(OUT_PNG, dpi=DPI, bbox_inches=None)  # bbox_inches=None to preserve exact pixel size
plt.savefig(OUT_SVG, bbox_inches="tight")

print(f"Saved: {OUT_PNG}")
print(f"Saved: {OUT_SVG}")
