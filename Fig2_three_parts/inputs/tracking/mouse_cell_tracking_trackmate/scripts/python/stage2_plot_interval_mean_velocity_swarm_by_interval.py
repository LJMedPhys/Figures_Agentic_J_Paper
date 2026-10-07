# Stage 2 (Plotting only): Swarmplot of per-interval mean velocities (one dot per interval)
# Input:  /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv
# Output: /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocity_swarm_by_interval.png/.svg

# NOTE: Do NOT re-import pandas/numpy/matplotlib/seaborn/scipy.stats/os; they are preloaded as pd/np/plt/sns/stats/os.

# -----------------------------
# Parameters
# -----------------------------
INPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv"
PROJECT_DIR = r"/app/data/projects/mouse_cell_tracking_trackmate"
FIG_DIR = os.path.join(PROJECT_DIR, "figures")
BASENAME = "interval_mean_velocity_swarm_by_interval"

Y_COL = "mean_velocity_um_per_min"
X_COL = "interval_label"

# Colorblind-safe blue/orange guidance: use blue for points
POINT_COLOR = "#0173B2"  # colorblind-safe blue
MEAN_LINE_COLOR = "black"

# -----------------------------
# Load
# -----------------------------
df = pd.read_csv(INPUT_CSV)

required = {X_COL, Y_COL}
missing = required - set(df.columns)
if missing:
    raise ValueError(
        f"Missing required columns in {INPUT_CSV}: {sorted(missing)}. Found: {list(df.columns)}"
    )

# Ensure ordering by FRAME if present; otherwise keep file order
if "FRAME" in df.columns:
    df = df.sort_values("FRAME").reset_index(drop=True)

# Coerce numeric
_df_y = pd.to_numeric(df[Y_COL], errors="coerce")
if _df_y.isna().any():
    bad_n = int(_df_y.isna().sum())
    raise ValueError(f"Found {bad_n} non-numeric/NaN values in {Y_COL} after coercion.")
df[Y_COL] = _df_y

# -----------------------------
# Summary + outliers (IQR) reporting (no removal)
# -----------------------------
overall_mean = float(df[Y_COL].mean())

q1 = float(df[Y_COL].quantile(0.25))
q3 = float(df[Y_COL].quantile(0.75))
iqr = q3 - q1
lower = q1 - 1.5 * iqr
upper = q3 + 1.5 * iqr
outlier_mask = (df[Y_COL] < lower) | (df[Y_COL] > upper)
outlier_count = int(outlier_mask.sum())

print(f"Loaded intervals: N={len(df)}")
print(f"Overall mean of {Y_COL}: {overall_mean:.6g} µm/min")
print(
    "IQR outlier check (threshold=1.5): "
    f"Q1={q1:.6g}, Q3={q3:.6g}, IQR={iqr:.6g}, lower={lower:.6g}, upper={upper:.6g}, outliers={outlier_count}"
)

# -----------------------------
# Plot
# -----------------------------
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_style("white")
sns.set_palette("colorblind")

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "axes.linewidth": 1.0,
})

fig, ax = plt.subplots(figsize=(10, 6))

# Swarmplot: one dot per interval; x is categorical interval label
sns.swarmplot(
    data=df,
    x=X_COL,
    y=Y_COL,
    color=POINT_COLOR,
    size=6,
    ax=ax,
)

# Overall mean line across full axis width
ax.axhline(overall_mean, color=MEAN_LINE_COLOR, linewidth=2.0, zorder=0)

ax.set_xlabel("Frame interval (t→t+1)", fontsize=14)
ax.set_ylabel("Mean velocity (µm/min)", fontsize=14)
ax.set_title("Per-interval mean velocity (TrackMate)", fontsize=16)

# Improve readability of many intervals
# NOTE: tick_params does NOT accept 'ha'. Use set_xticklabels (or plt.setp) instead.
ax.tick_params(axis="x", rotation=45)
for lbl in ax.get_xticklabels():
    lbl.set_horizontalalignment("right")

# Clean spines for publication
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Optional: annotate mean value
ax.text(
    0.99,
    overall_mean,
    f" mean={overall_mean:.2f}",
    transform=ax.get_yaxis_transform(),
    ha="right",
    va="bottom",
    fontsize=12,
    color=MEAN_LINE_COLOR,
)

png_path = os.path.join(FIG_DIR, f"{BASENAME}.png")
svg_path = os.path.join(FIG_DIR, f"{BASENAME}.svg")

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(svg_path, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")
