# Stage 2 (Plotting only)
# Narrow compact single-category beeswarm plot of per-interval mean velocities.
# Reads:  /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv
# Writes: /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_beeswarm_narrow_y0_tick1_meansegment.png
#         /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_beeswarm_narrow_y0_tick1_meansegment.svg

# NOTE: Do not re-import pandas/numpy/matplotlib/seaborn/scipy.stats/os; they are pre-initialized.

# ----------------------------
# Config
# ----------------------------
INPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv"
PROJECT_DIR = r"/app/data/projects/mouse_cell_tracking_trackmate"
FIG_DIR = os.path.join(PROJECT_DIR, "figures")

OUT_PNG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_narrow_y0_tick1_meansegment.png")
OUT_SVG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_narrow_y0_tick1_meansegment.svg")

Y_COL = "mean_velocity_um_per_min"

# Narrow figure: clearly narrower than prior compact plots
FIG_W, FIG_H = 1.45, 4.2

POINT_COLOR = "#0173B2"  # colorblind-safe blue

# ----------------------------
# Load
# ----------------------------
df = pd.read_csv(INPUT_CSV)

if Y_COL not in df.columns:
    raise ValueError(f"Expected column '{Y_COL}' not found. Columns: {list(df.columns)}")

# Keep finite values only
vals = pd.to_numeric(df[Y_COL], errors="coerce")
df = df.loc[vals.notna() & np.isfinite(vals)].copy()

if df.empty:
    raise ValueError("No valid numeric values found for plotting after filtering NaN/inf.")

# Single category for swarm packing
cat_label = ""
df["_category"] = cat_label

# Overall mean across per-interval means (not recomputing any hypothesis tests)
overall_mean = float(df[Y_COL].mean())
print(f"Overall mean of {Y_COL} across intervals: {overall_mean:.6f} µm/min (n={len(df)})")

# ----------------------------
# Plot
# ----------------------------
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_style("white")
sns.set_palette("colorblind")

plt.figure(figsize=(FIG_W, FIG_H))
ax = sns.swarmplot(
    data=df,
    x="_category",
    y=Y_COL,
    size=5,
    color=POINT_COLOR,
    edgecolor=None,
    linewidth=0,
)

# Axis labels/ticks: single central axis, no x ticks/labels
ax.set_xlabel("")
ax.set_xticks([])
ax.set_xticklabels([])

ax.set_ylabel("Mean velocity (µm/min)", fontsize=14)
ax.tick_params(axis="y", labelsize=12)

# y-axis starts at 0 with major ticks every 1.0
ymax_data = float(df[Y_COL].max())
ymax = max(1.0, np.ceil(ymax_data))
ax.set_ylim(0, ymax)
ax.set_yticks(np.arange(0, ymax + 1e-9, 1.0))

# Minimal styling: no grid; left/bottom spines only
ax.grid(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1.0)
ax.spines["bottom"].set_linewidth(1.0)

# Short black overall-mean segment spanning just the dot column region
# For a single category at x=0, draw a short segment centered at 0.
segment_half_width = 0.12
ax.plot(
    [-segment_half_width, segment_half_width],
    [overall_mean, overall_mean],
    color="black",
    linewidth=2.0,
    solid_capstyle="butt",
    zorder=5,
)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.savefig(OUT_SVG, bbox_inches="tight")
plt.close()

print(f"Saved: {OUT_PNG}")
print(f"Saved: {OUT_SVG}")
