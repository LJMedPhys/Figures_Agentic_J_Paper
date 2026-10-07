# Stage 2 (plotting-only): Beeswarm plot of per-interval mean velocities with short overall-mean segment
# Inputs:
#   - /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv
#   - /app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv (preferred for overall mean)
# Outputs:
#   - /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_step_velocity_beeswarm_meansegment_y0_tick1.png
#   - /app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_step_velocity_beeswarm_meansegment_y0_tick1.svg
# Notes:
#   - No statistical testing performed here (plotting-only).
#   - Overall mean is computed from all step velocities if available; otherwise from interval means.

# Environment presets: pd, np, plt, sns, stats, os are already imported.

project_dir = r"/app/data/projects/mouse_cell_tracking_trackmate"
input_interval_csv = os.path.join(project_dir, "figures", "interval_mean_velocities.csv")
input_stats_csv = os.path.join(project_dir, "data", "Statistics_Results.csv")
fig_dir = os.path.join(project_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, "interval_mean_step_velocity_beeswarm_meansegment_y0_tick1.png")
out_svg = os.path.join(fig_dir, "interval_mean_step_velocity_beeswarm_meansegment_y0_tick1.svg")

# -------------------------
# Load interval means
# -------------------------
df_int = pd.read_csv(input_interval_csv)
required_col = "mean_velocity_um_per_min"
if required_col not in df_int.columns:
    raise ValueError(f"Expected column '{required_col}' in {input_interval_csv}. Found: {list(df_int.columns)}")

interval_means = df_int[required_col].astype(float).dropna().values
n_intervals = len(interval_means)
if n_intervals == 0:
    raise ValueError("No interval mean velocities found after dropping NaNs.")

# -------------------------
# Determine overall mean (prefer step velocities)
# -------------------------
mean_source = "interval_means"
overall_mean = float(np.mean(interval_means))

if os.path.exists(input_stats_csv):
    try:
        df_steps = pd.read_csv(input_stats_csv)
        step_col = "STEP_VELOCITY_UM_PER_MIN"
        if step_col in df_steps.columns:
            step_vals = df_steps[step_col].astype(float).dropna().values
            if len(step_vals) > 0:
                overall_mean = float(np.mean(step_vals))
                mean_source = "all_step_velocities"
    except Exception as e:
        print(f"WARNING: Could not parse {input_stats_csv} for overall mean; falling back to interval means. Reason: {e}")

print(f"N intervals plotted: {n_intervals}")
print(f"Overall mean velocity used (um/min): {overall_mean:.6f} (source={mean_source})")

# -------------------------
# Plot styling (publication)
# -------------------------
sns.set_theme(style="white")
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

# Prepare plotting dataframe: single category
plot_df = pd.DataFrame({
    "category": ["Interval mean"] * n_intervals,
    "velocity_um_per_min": interval_means
})

fig, ax = plt.subplots(figsize=(4.0, 6.0))

# Beeswarm (swarmplot) with blue points
blue = "#0173B2"
# Use a moderate size to reduce overlap; swarmplot will pack points
sns.swarmplot(
    data=plot_df,
    x="category",
    y="velocity_um_per_min",
    ax=ax,
    color=blue,
    size=6,
    edgecolor="none",
    linewidth=0,
)

# Remove x-axis ticks/labels
ax.set_xlabel("")
ax.set_xticks([])
ax.set_xticklabels([])

# Axis labels
ax.set_ylabel("Mean step velocity per interval (µm/min)")

# Y-axis starts at 0 with major ticks every 1.0
# Determine upper limit with a small headroom
ymax_data = float(np.nanmax(interval_means))
ymax = max(ymax_data, overall_mean)
headroom = 0.08 * (ymax if ymax > 0 else 1.0)
ymax_lim = ymax + headroom
ax.set_ylim(0, ymax_lim)

# Major ticks every 1.0
from matplotlib.ticker import MultipleLocator
ax.yaxis.set_major_locator(MultipleLocator(1.0))

# No gridlines
ax.grid(False)

# Only left + bottom spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# -------------------------
# Short overall-mean segment spanning swarm width
# -------------------------
# Swarm points are drawn at x positions around 0 (single category). We can infer the swarm width
# from the PathCollections in the axis.

x_min, x_max = None, None
for coll in ax.collections:
    offsets = getattr(coll, "get_offsets", lambda: None)()
    if offsets is None or len(offsets) == 0:
        continue
    xs = np.asarray(offsets)[:, 0]
    if xs.size == 0:
        continue
    x_min = float(np.min(xs)) if x_min is None else min(x_min, float(np.min(xs)))
    x_max = float(np.max(xs)) if x_max is None else max(x_max, float(np.max(xs)))

# Fallback if offsets not found
if x_min is None or x_max is None:
    x_min, x_max = -0.2, 0.2

# Add a small padding so the segment slightly exceeds the swarm width
pad = 0.03
seg_x0 = x_min - pad
seg_x1 = x_max + pad

ax.plot([seg_x0, seg_x1], [overall_mean, overall_mean], color="black", linewidth=2.0, solid_capstyle="butt")

# Tight layout and save
plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
