# Plotting-only script (Stage 2)
# Creates a compact single-category jitter/strip plot of per-interval mean velocities.
# Prefers precomputed interval means CSV; otherwise computes from per-step Statistics_Results.csv.

# Environment presets assumed: pd, np, plt, sns, stats, os

project_dir = r"/app/data/projects/mouse_cell_tracking_trackmate"
fig_dir = os.path.join(project_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

interval_means_path = os.path.join(fig_dir, "interval_mean_velocities.csv")
steps_path = os.path.join(project_dir, "data", "Statistics_Results.csv")

out_png = os.path.join(fig_dir, "interval_mean_step_velocity_compact_jitter.png")
out_svg = os.path.join(fig_dir, "interval_mean_step_velocity_compact_jitter.svg")

# -------------------------
# Load interval means (preferred) or compute
# -------------------------
interval_df = None
source_used = None

if os.path.exists(interval_means_path):
    interval_df = pd.read_csv(interval_means_path)
    source_used = "interval_mean_velocities.csv"
else:
    # Fallback: compute per-interval means from per-step velocities
    # Robust read in case of occasional malformed lines
    try:
        steps_df = pd.read_csv(steps_path)
    except Exception:
        steps_df = pd.read_csv(steps_path, engine="python", on_bad_lines="skip")

    required = {"FRAME", "FRAME_NEXT", "STEP_VELOCITY_UM_PER_MIN"}
    missing = required - set(steps_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {steps_path}: {sorted(missing)}. Found: {list(steps_df.columns)}")

    interval_df = (
        steps_df
        .dropna(subset=["FRAME", "FRAME_NEXT", "STEP_VELOCITY_UM_PER_MIN"])
        .groupby(["FRAME", "FRAME_NEXT"], as_index=False)
        .agg(
            mean_velocity_um_per_min=("STEP_VELOCITY_UM_PER_MIN", "mean"),
            n_steps=("STEP_VELOCITY_UM_PER_MIN", "size"),
            sd_velocity_um_per_min=("STEP_VELOCITY_UM_PER_MIN", "std"),
        )
    )
    interval_df["interval_label"] = interval_df["FRAME"].astype(int).astype(str) + "→" + interval_df["FRAME_NEXT"].astype(int).astype(str)
    source_used = "computed_from_Statistics_Results.csv"

# Validate expected column
if "mean_velocity_um_per_min" not in interval_df.columns:
    raise ValueError(
        "Expected column 'mean_velocity_um_per_min' in interval means table. "
        f"Found columns: {list(interval_df.columns)}"
    )

# -------------------------
# Compute overall mean across ALL steps (preferred for consistency)
# -------------------------
try:
    steps_df_for_mean = pd.read_csv(steps_path)
except Exception:
    steps_df_for_mean = pd.read_csv(steps_path, engine="python", on_bad_lines="skip")

if "STEP_VELOCITY_UM_PER_MIN" not in steps_df_for_mean.columns:
    raise ValueError(
        "Expected column 'STEP_VELOCITY_UM_PER_MIN' in Statistics_Results.csv. "
        f"Found columns: {list(steps_df_for_mean.columns)}"
    )

step_vel = steps_df_for_mean["STEP_VELOCITY_UM_PER_MIN"].dropna().astype(float)
overall_mean_steps = float(step_vel.mean())

# IQR outlier reporting (no removal)
q1, q3 = np.percentile(step_vel, [25, 75])
iqr = q3 - q1
low = q1 - 1.5 * iqr
high = q3 + 1.5 * iqr
n_outliers = int(((step_vel < low) | (step_vel > high)).sum())

print(f"Interval means source used: {source_used}")
print(f"N intervals plotted: {interval_df.shape[0]}")
print(f"Overall mean velocity across ALL steps (µm/min): {overall_mean_steps:.6f}")
print(f"IQR outliers in step velocities (threshold=1.5): {n_outliers} / {step_vel.shape[0]}")

# -------------------------
# Plot: compact single-category jitter
# -------------------------
sns.set_style("white")
sns.set_palette("colorblind")

plt.figure(figsize=(4.2, 6.0))
ax = plt.gca()

# Single x category at 0; jitter points around it
x0 = np.zeros(interval_df.shape[0], dtype=float)
# Deterministic jitter for reproducibility
rng = np.random.default_rng(0)
jitter = rng.uniform(-0.08, 0.08, size=interval_df.shape[0])

ax.scatter(
    x0 + jitter,
    interval_df["mean_velocity_um_per_min"].astype(float).values,
    s=35,
    color="#0173B2",  # colorblind-safe blue
    edgecolor="black",
    linewidth=0.4,
    alpha=0.9,
    zorder=3,
)

# Overall mean line (across all steps)
ax.axhline(overall_mean_steps, color="black", linewidth=2.0, zorder=2)

# Formatting: no x labels/title, no grid
ax.set_xlim(-0.35, 0.35)
ax.set_xticks([])
ax.set_xlabel("")
ax.set_ylabel("Velocity (µm/min)", fontsize=14)
ax.tick_params(axis="y", labelsize=12)
ax.grid(False)

# Spine styling
for spine in ax.spines.values():
    spine.set_linewidth(1.0)

ax.set_title(
    "Per-interval mean step velocity\n(TrackMate LAP tracking; mean line = all steps)",
    fontsize=16,
    pad=10,
)

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
plt.close()

print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
