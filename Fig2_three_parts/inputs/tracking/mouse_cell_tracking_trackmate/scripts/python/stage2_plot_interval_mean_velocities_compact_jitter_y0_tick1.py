# Stage 2 (plotting-only)
# Compact single-axis jitter plot of per-interval mean velocities with publication styling.

# Environment presets: pd, np, plt, sns, stats, os are already imported.

project_dir = r"/app/data/projects/mouse_cell_tracking_trackmate"
input_csv = r"/app/data/projects/mouse_cell_tracking_trackmate/figures/interval_mean_velocities.csv"
fig_dir = os.path.join(project_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, "interval_mean_step_velocity_compact_jitter_y0_tick1.png")
out_svg = os.path.join(fig_dir, "interval_mean_step_velocity_compact_jitter_y0_tick1.svg")

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(input_csv)

required_cols = ["mean_velocity_um_per_min"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in {input_csv}: {missing}. Found columns: {list(df.columns)}")

vals = pd.to_numeric(df["mean_velocity_um_per_min"], errors="coerce").dropna()
if vals.empty:
    raise ValueError("No valid numeric values found in mean_velocity_um_per_min.")

# Overall mean across intervals (points are interval means)
overall_mean = float(vals.mean())
print(f"N intervals = {len(vals)}")
print(f"Overall mean of interval mean velocities (um/min) = {overall_mean:.6g}")

# -----------------------------
# Plot styling
# -----------------------------
sns.set_theme(style="ticks")
sns.set_palette("colorblind")

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.0,
})

# Single central x with jitter
rng = np.random.default_rng(0)
x0 = 0.0
jitter = rng.uniform(-0.08, 0.08, size=len(vals))
x = x0 + jitter

fig, ax = plt.subplots(figsize=(4.0, 6.0))

# Points (blue)
ax.scatter(
    x,
    vals.values,
    s=35,
    color="#0173B2",  # colorblind-safe blue
    edgecolor="none",
    alpha=0.9,
    zorder=3,
)

# Mean line (black)
ax.axhline(overall_mean, color="black", linewidth=2.0, zorder=2)

# Axis labels and limits
ax.set_ylabel("Mean step velocity per interval (µm/min)")
ax.set_xlabel("")
ax.set_xlim(-0.35, 0.35)

# y-axis starts at 0
upper = float(np.nanmax(vals.values))
# Add headroom so top point/mean line isn't cramped
ax.set_ylim(0, upper * 1.08 if upper > 0 else 1)

# y major ticks increment by 1.0
from matplotlib.ticker import MultipleLocator
ax.yaxis.set_major_locator(MultipleLocator(1.0))

# Remove x ticks/labels entirely
ax.set_xticks([])
ax.set_xticklabels([])

# Remove all gridlines
ax.grid(False)

# Show only left/bottom spines (axis lines)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# Ensure ticks only on left/bottom
ax.tick_params(axis="both", which="both", direction="out", length=4, width=1.0, top=False, right=False)

# Save
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
