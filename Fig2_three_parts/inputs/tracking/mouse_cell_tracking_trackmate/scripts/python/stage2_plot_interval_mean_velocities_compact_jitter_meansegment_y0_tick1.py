# Stage 2 (plotting-only)
# Interval mean velocities: compact single-axis jitter (strip) with short overall-mean segment.
# Reads per-interval mean velocities from figures/interval_mean_velocities.csv.
# Uses overall mean from per-step velocities in data/Statistics_Results.csv if available; otherwise mean of interval means.

# NOTE: Do NOT re-import pandas/numpy/matplotlib/seaborn/scipy.stats/os; they are pre-initialized in the environment.

# ---------------------------
# Paths
# ---------------------------
project_dir = r"/app/data/projects/mouse_cell_tracking_trackmate"
points_csv = os.path.join(project_dir, "figures", "interval_mean_velocities.csv")
step_stats_csv = os.path.join(project_dir, "data", "Statistics_Results.csv")

fig_dir = os.path.join(project_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, "interval_mean_step_velocity_compact_jitter_meansegment_y0_tick1.png")
out_svg = os.path.join(fig_dir, "interval_mean_step_velocity_compact_jitter_meansegment_y0_tick1.svg")

# ---------------------------
# Load data
# ---------------------------
df = pd.read_csv(points_csv)

required_cols = ["mean_velocity_um_per_min"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in {points_csv}: {missing}. Found columns: {list(df.columns)}")

vel = df["mean_velocity_um_per_min"].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
if vel.empty:
    raise ValueError("No valid mean_velocity_um_per_min values found after cleaning.")

# ---------------------------
# Determine overall mean (prefer per-step velocities)
# ---------------------------
mean_source = "interval_means"
overall_mean = float(vel.mean())

if os.path.exists(step_stats_csv):
    try:
        # Robust read: tolerate occasional malformed lines from prior pipeline stages
        step_df = pd.read_csv(step_stats_csv, on_bad_lines='skip', engine='python')
        if "STEP_VELOCITY_UM_PER_MIN" in step_df.columns:
            step_vel = step_df["STEP_VELOCITY_UM_PER_MIN"].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
            if not step_vel.empty:
                overall_mean = float(step_vel.mean())
                mean_source = "all_steps"
    except Exception as e:
        print(f"WARNING: Could not parse {step_stats_csv}; falling back to interval mean-of-means. Reason: {e}")

print(f"Overall mean velocity used ({mean_source}) = {overall_mean:.6f} µm/min")

# ---------------------------
# Outlier reporting (IQR) - report only, do not remove
# ---------------------------
q1, q3 = np.percentile(vel.values, [25, 75])
iqr = q3 - q1
lower = q1 - 1.5 * iqr
upper = q3 + 1.5 * iqr
n_outliers = int(((vel < lower) | (vel > upper)).sum())
print(f"IQR outliers (threshold=1.5): {n_outliers} of {len(vel)} interval means")

# ---------------------------
# Plot styling
# ---------------------------
sns.set_theme(style="white")
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

# Single central x with horizontal jitter
x0 = 0
x = np.full(len(vel), x0, dtype=float)

# Deterministic jitter for reproducibility
rng = np.random.default_rng(7)
jitter_width = 0.12
x_j = x + rng.uniform(-jitter_width, jitter_width, size=len(vel))

# Colors (colorblind-safe): blue for points
point_color = "#0173B2"

fig, ax = plt.subplots(figsize=(4.0, 6.0))

# Scatter points
ax.scatter(
    x_j,
    vel.values,
    s=28,
    color=point_color,
    edgecolor="none",
    alpha=0.9,
    zorder=3,
)

# Short mean segment spanning the point column region (not full width)
segment_halfwidth = jitter_width * 1.15
ax.hlines(
    overall_mean,
    x0 - segment_halfwidth,
    x0 + segment_halfwidth,
    colors="black",
    linewidth=2.0,
    zorder=4,
)

# Axes formatting
ax.set_xlim(x0 - 0.35, x0 + 0.35)
ax.set_ylim(bottom=0)

# y ticks every 1.0
ymax = float(np.nanmax(vel.values))
# Ensure at least one tick above max for headroom
upper_tick = np.ceil(max(ymax, overall_mean) + 0.25)
ax.set_yticks(np.arange(0, upper_tick + 1e-9, 1.0))

ax.set_ylabel("Mean step velocity per interval (µm/min)")

# No x-axis ticks/labels
ax.set_xticks([])
ax.set_xlabel("")

# No gridlines
ax.grid(False)

# Only left/bottom spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# Tight layout and save
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
