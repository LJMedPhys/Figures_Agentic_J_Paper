# Stage 2 (Plotting only)
# Publication-quality single-category beeswarm plot of per-interval mean velocities.
# Reads: /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv
# Writes: figures/interval_mean_velocity_beeswarm_padded_ylabel22_ytickfs18_ylim3_meansegment.(png|svg)

# NOTE: Do not re-import pandas/numpy/matplotlib/seaborn/scipy.stats/os; they are preloaded in the environment.

# -----------------------------
# Parameters (styling + outputs)
# -----------------------------
PROJECT_DIR = r"/app/data/projects/mouse_cell_tracking_trackmate"
INPUT_CSV = os.path.join(PROJECT_DIR, "figures", "interval_mean_velocities.csv")
FIG_DIR = os.path.join(PROJECT_DIR, "figures")

OUT_PNG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_padded_ylabel22_ytickfs18_ylim3_meansegment.png")
OUT_SVG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_padded_ylabel22_ytickfs18_ylim3_meansegment.svg")

Y_COL = "mean_velocity_um_per_min"
YLABEL = "Mean velocity (µm/min)"

# Match prior reference styling as closely as possible
DOT_COLOR = "#0173B2"  # colorblind-safe blue
DOT_SIZE = 5           # keep same as reference plot
MEAN_LINE_COLOR = "black"
MEAN_LINE_WIDTH = 2.0

# Figure size: match prior 530x600 px at 300 dpi
DPI = 300
FIG_W_IN = 530 / DPI
FIG_H_IN = 600 / DPI

# Axis limits/ticks requested
YMIN, YMAX = 0, 3
YTICKS = [0, 1, 2, 3]
YTICK_FONTSIZE = 18
YLABEL_FONTSIZE = 22

# Mean segment geometry (short segment spanning only dot column)
SEG_X0, SEG_X1 = -0.08, 0.08

# Padding/margins (tuned to resemble prior padded plot)
LEFT_PAD = 0.28
RIGHT_PAD = 0.95
BOTTOM_PAD = 0.08
TOP_PAD = 0.98

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(INPUT_CSV)

if Y_COL not in df.columns:
    raise ValueError(f"Expected column '{Y_COL}' not found. Columns: {list(df.columns)}")

# Drop missing values defensively
vals = df[Y_COL].astype(float).dropna().values
if vals.size == 0:
    raise ValueError("No valid velocity values found after dropping NaNs.")

overall_mean = float(np.mean(vals))
print(f"Overall mean of per-interval mean velocities ({Y_COL}): {overall_mean:.6f} µm/min")

# -----------------------------
# Plot
# -----------------------------
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_style("white")
sns.set_palette("colorblind")

plt.figure(figsize=(FIG_W_IN, FIG_H_IN))
ax = plt.gca()

# Single-category swarm: use constant x category
plot_df = df[[Y_COL]].copy()
plot_df["category"] = ""

sns.swarmplot(
    data=plot_df,
    x="category",
    y=Y_COL,
    color=DOT_COLOR,
    size=DOT_SIZE,
    ax=ax
)

# Overall mean short segment centered on dot column
ax.plot([SEG_X0, SEG_X1], [overall_mean, overall_mean], color=MEAN_LINE_COLOR, linewidth=MEAN_LINE_WIDTH)

# Axes formatting
ax.set_xlabel("")
ax.set_ylabel(YLABEL, fontsize=YLABEL_FONTSIZE)

ax.set_ylim(YMIN, YMAX)
ax.set_yticks(YTICKS)

# Remove x ticks/labels entirely
ax.set_xticks([])
ax.set_xticklabels([])

# Minimal spines: left + bottom only
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

# No grid
ax.grid(False)

# Y tick labels + visible tick marks on left, outward
ax.tick_params(
    axis="y",
    which="both",
    left=True,
    right=False,
    direction="out",
    length=6,
    width=1.2,
    labelsize=YTICK_FONTSIZE
)

# Ensure x tick marks are off
ax.tick_params(axis="x", bottom=False, top=False, labelbottom=False)

# Tight layout with controlled padding to match reference
plt.subplots_adjust(left=LEFT_PAD, right=RIGHT_PAD, bottom=BOTTOM_PAD, top=TOP_PAD)

plt.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
plt.savefig(OUT_SVG, bbox_inches="tight")
plt.close()

print(f"Saved: {OUT_PNG}")
print(f"Saved: {OUT_SVG}")
