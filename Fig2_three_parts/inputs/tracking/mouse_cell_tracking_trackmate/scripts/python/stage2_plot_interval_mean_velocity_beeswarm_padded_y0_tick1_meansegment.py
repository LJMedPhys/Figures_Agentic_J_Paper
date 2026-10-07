# Stage 2 (Plotting only)
# Publication-quality beeswarm (single-category swarmplot) of per-interval mean velocities.
# Reads:  /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv
# Writes: /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_beeswarm_padded_y0_tick1_meansegment.png
#         /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_beeswarm_padded_y0_tick1_meansegment.svg

# NOTE: No statistical tests are performed here (Stage 2 rule). Overall mean is descriptive.

import os

# ----------------------------
# Paths
# ----------------------------
PROJECT_DIR = r"/app/data/projects/mouse_cell_tracking_trackmate"
INPUT_CSV = os.path.join(PROJECT_DIR, "figures", "interval_mean_velocities.csv")
FIG_DIR = os.path.join(PROJECT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

OUT_PNG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_padded_y0_tick1_meansegment.png")
OUT_SVG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_padded_y0_tick1_meansegment.svg")

# ----------------------------
# Load data
# ----------------------------
df = pd.read_csv(INPUT_CSV)

required_cols = ["mean_velocity_um_per_min"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}. Found columns: {list(df.columns)}")

# Clean numeric column
vel = pd.to_numeric(df["mean_velocity_um_per_min"], errors="coerce")
vel = vel.replace([np.inf, -np.inf], np.nan).dropna()

if vel.empty:
    raise ValueError("No valid numeric values found in mean_velocity_um_per_min after cleaning.")

overall_mean = float(vel.mean())

# ----------------------------
# Outlier reporting (IQR) - report only, do not remove
# ----------------------------
q1 = float(vel.quantile(0.25))
q3 = float(vel.quantile(0.75))
iqr = q3 - q1
if iqr == 0:
    outlier_mask = pd.Series([False] * len(vel), index=vel.index)
else:
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    outlier_mask = (vel < lower) | (vel > upper)

n_outliers = int(outlier_mask.sum())
print(f"N intervals (points) = {len(vel)}")
print(f"Overall mean(mean_velocity_um_per_min) = {overall_mean:.6f} µm/min")
print(f"IQR outliers (threshold=1.5): {n_outliers} / {len(vel)}")

# ----------------------------
# Plot styling
# ----------------------------
sns.set_style("white")
sns.set_palette("colorblind")

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "axes.linewidth": 1.0,
})

# Slightly wider than prior narrow version + explicit left/right padding.
# (Exact pixel padding depends on backend; we enforce via subplots_adjust.)
fig, ax = plt.subplots(figsize=(3.0, 6.0))

plot_df = pd.DataFrame({
    "group": [""] * len(vel),
    "mean_velocity_um_per_min": vel.values
})

# Beeswarm (single category)
sns.swarmplot(
    data=plot_df,
    x="group",
    y="mean_velocity_um_per_min",
    ax=ax,
    size=5,
    color="#0173B2",  # colorblind-safe blue
    edgecolor=None,
    linewidth=0,
)

# Overall mean as short segment centered on the dot column
# For a single category, x is at 0. Draw a short segment around it.
segment_half_width = 0.12
ax.hlines(
    overall_mean,
    xmin=0 - segment_half_width,
    xmax=0 + segment_half_width,
    colors="black",
    linewidth=2.0,
    zorder=5,
)

# Axes formatting
ax.set_xlabel("")
ax.set_ylabel("Mean velocity (µm/min)", fontsize=24)
ax.set_title("")

# y-axis from 0 with 1.0 tick spacing
ax.set_ylim(bottom=0)
ymax = float(np.nanmax(vel.values))
# Ensure at least one tick above max
upper_tick = np.ceil(max(ymax, overall_mean) / 1.0) * 1.0
ax.set_yticks(np.arange(0, upper_tick + 1.0, 1.0))

# Remove x ticks/labels
ax.set_xticks([])
ax.set_xticklabels([])

# Minimal spines: left/bottom only
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1.0)
ax.spines["bottom"].set_linewidth(1.0)

# No grid
ax.grid(False)

# Add explicit left/right padding so points don't collide with axes.
# Also widen x-limits slightly to create whitespace around the swarm.
ax.set_xlim(-0.6, 0.6)

# Subplot padding: tuned to create visible whitespace on both sides.
# (Approx ~50 px depends on final rendered size; this is a robust visual padding.)
fig.subplots_adjust(left=0.28, right=0.92, bottom=0.08, top=0.98)

# Save
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_PNG}")
print(f"Saved: {OUT_SVG}")
