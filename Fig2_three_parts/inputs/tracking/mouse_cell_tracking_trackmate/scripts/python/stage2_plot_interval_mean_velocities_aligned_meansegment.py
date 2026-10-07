# Stage 2 (plotting-only)
# Aligned (no-jitter) interval-mean velocity plot with short overall-mean segment.
# Inputs:
#   - figures/interval_mean_velocities.csv (one row per interval)
#   - data/Statistics_Results.csv (optional; used to compute overall mean across all steps)
# Outputs:
#   - figures/interval_mean_step_velocity_aligned_meansegment.png (300 DPI)
#   - figures/interval_mean_step_velocity_aligned_meansegment.svg

import os

# ---------------------------
# Paths
# ---------------------------
project_dir = r"/app/data/projects/mouse_cell_tracking_trackmate"
interval_csv = os.path.join(project_dir, "figures", "interval_mean_velocities.csv")
steps_stats_csv = os.path.join(project_dir, "data", "Statistics_Results.csv")
fig_dir = os.path.join(project_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, "interval_mean_step_velocity_aligned_meansegment.png")
out_svg = os.path.join(fig_dir, "interval_mean_step_velocity_aligned_meansegment.svg")

# ---------------------------
# Load interval means
# ---------------------------
df_int = pd.read_csv(interval_csv)

required_cols = {"mean_velocity_um_per_min"}
missing = required_cols - set(df_int.columns)
if missing:
    raise ValueError(f"interval_mean_velocities.csv missing required columns: {sorted(missing)}. Found: {list(df_int.columns)}")

y = pd.to_numeric(df_int["mean_velocity_um_per_min"], errors="coerce").dropna()
if y.empty:
    raise ValueError("No valid numeric values found in mean_velocity_um_per_min.")

# ---------------------------
# Determine overall mean velocity
# Prefer overall mean across all steps from Statistics_Results.csv if available.
# ---------------------------
mean_source = "interval_means"
overall_mean = float(y.mean())

if os.path.exists(steps_stats_csv):
    try:
        # Robust parsing in case file has inconsistent columns (prior tokenizing issues).
        df_steps = pd.read_csv(steps_stats_csv, engine="python", on_bad_lines="skip")
        if "STEP_VELOCITY_UM_PER_MIN" in df_steps.columns:
            step_v = pd.to_numeric(df_steps["STEP_VELOCITY_UM_PER_MIN"], errors="coerce").dropna()
            if len(step_v) > 0:
                overall_mean = float(step_v.mean())
                mean_source = "all_steps"
    except Exception as e:
        print(f"WARNING: Could not parse Statistics_Results.csv for overall mean across steps; using mean of interval means instead. Reason: {e}")

print(f"Overall mean velocity used = {overall_mean:.6f} µm/min (source: {mean_source}).")

# ---------------------------
# Outlier reporting (IQR method) - report only, do not remove
# ---------------------------
q1 = np.percentile(y, 25)
q3 = np.percentile(y, 75)
iqr = q3 - q1
lower = q1 - 1.5 * iqr
upper = q3 + 1.5 * iqr
n_outliers = int(((y < lower) | (y > upper)).sum())
print(f"IQR outliers (threshold=1.5): {n_outliers} of {len(y)} interval means.")

# ---------------------------
# Plot styling (publication)
# ---------------------------
sns.set_style("white")
sns.set_palette("colorblind")
plt.rcParams.update({
    "font.size": 12,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})

fig, ax = plt.subplots(figsize=(4.0, 6.0))

# Single central x position (no jitter)
x0 = 0.0
x = np.full(len(y), x0)

# Blue points aligned vertically
ax.scatter(
    x, y,
    s=35,
    color="#0173B2",  # colorblind-safe blue
    edgecolor="none",
    alpha=0.9,
    zorder=3
)

# Short mean segment crossing only the point column region
segment_half_width = 0.12
ax.plot(
    [x0 - segment_half_width, x0 + segment_half_width],
    [overall_mean, overall_mean],
    color="black",
    linewidth=2.0,
    solid_capstyle="butt",
    zorder=4
)

# Axes labels and limits
ax.set_ylabel("Mean step velocity per interval (µm/min)", fontsize=14)
ax.set_xlabel("")

# No x-axis ticks/labels
ax.set_xticks([])
ax.set_xlim(-0.5, 0.5)

# y-axis starts at 0 with tick spacing 1.0
ymax = float(max(y.max(), overall_mean))
y_upper = max(1.0, np.ceil((ymax + 0.25) / 1.0) * 1.0)
ax.set_ylim(0, y_upper)
ax.set_yticks(np.arange(0, y_upper + 1e-9, 1.0))
ax.tick_params(axis="y", labelsize=12)

# No gridlines
ax.grid(False)

# Only left and bottom spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1.0)
ax.spines["bottom"].set_linewidth(1.0)

# Title (kept minimal; can be edited in SVG)
ax.set_title("Interval mean velocities", fontsize=16, pad=10)

plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
