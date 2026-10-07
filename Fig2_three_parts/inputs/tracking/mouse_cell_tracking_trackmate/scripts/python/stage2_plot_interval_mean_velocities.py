# Stage 2 (plotting-only): Interval mean step velocities from TrackMate steps
# Reads Statistics_Results.csv (per-step rows), aggregates mean velocity per frame interval,
# and generates publication-quality strip/jitter plot.

# Environment presets: pd, np, plt, sns, stats, os are already imported.

project_dir = r"/app/data/projects/mouse_cell_tracking_trackmate"
input_csv = os.path.join(project_dir, "data", "Statistics_Results.csv")
fig_dir = os.path.join(project_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, "interval_mean_step_velocity_jitter.png")
out_svg = os.path.join(fig_dir, "interval_mean_step_velocity_jitter.svg")
out_interval_csv = os.path.join(fig_dir, "interval_mean_velocities.csv")


def _read_stats_csv_robust(path):
    """Robustly read CSV even if prior pipeline appended non-tabular lines.

    Strategy:
      1) Try normal pd.read_csv.
      2) If it fails, fall back to python engine with on_bad_lines='skip'.
      3) If still problematic, manually filter lines that look like data rows.
    """
    try:
        return pd.read_csv(path)
    except Exception as e1:
        print(f"Standard read_csv failed: {e1}")

    try:
        return pd.read_csv(path, engine="python", on_bad_lines="skip")
    except Exception as e2:
        print(f"Fallback read_csv(engine='python', on_bad_lines='skip') failed: {e2}")

    # Manual filter: keep header + lines with consistent comma count
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()

    if not lines:
        raise ValueError("Empty CSV file")

    header = lines[0]
    n_commas = header.count(",")
    kept = [header]
    for ln in lines[1:]:
        if ln.count(",") == n_commas:
            kept.append(ln)

    from io import StringIO
    return pd.read_csv(StringIO("\n".join(kept)))


# -------------------- Load & validate --------------------
df = _read_stats_csv_robust(input_csv)
print("Loaded columns:", list(df.columns))
print("Head:\n", df.head())

# Identify required columns with allowed aliases
col_from = None
col_to = None
col_vel = None

for c in ["from_FRAME", "FRAME", "FROM_FRAME", "From_FRAME"]:
    if c in df.columns:
        col_from = c
        break

for c in ["to_FRAME", "FRAME_NEXT", "TO_FRAME", "Frame_Next", "FRAME_TO"]:
    if c in df.columns:
        col_to = c
        break

for c in ["STEP_VELOCITY_UM_PER_MIN", "velocity_um_per_min", "STEP_VELOCITY", "VELOCITY_UM_PER_MIN"]:
    if c in df.columns:
        col_vel = c
        break

if col_from is None or col_to is None or col_vel is None:
    raise KeyError(
        "Could not find required columns. "
        f"Found from={col_from}, to={col_to}, vel={col_vel}. "
        f"Available columns: {list(df.columns)}"
    )

# Coerce numeric and drop invalid rows
work = df[[col_from, col_to, col_vel]].copy()
work[col_from] = pd.to_numeric(work[col_from], errors="coerce")
work[col_to] = pd.to_numeric(work[col_to], errors="coerce")
work[col_vel] = pd.to_numeric(work[col_vel], errors="coerce")

n0 = len(work)
work = work.dropna(subset=[col_from, col_to, col_vel])
print(f"Rows kept after dropping NA in interval/velocity: {len(work)}/{n0}")

# Ensure integer frames
work[col_from] = work[col_from].astype(int)
work[col_to] = work[col_to].astype(int)

# -------------------- Aggregate per interval --------------------
interval_stats = (
    work.groupby([col_from, col_to], as_index=False)
        .agg(
            mean_velocity_um_per_min=(col_vel, "mean"),
            n_steps=(col_vel, "size"),
            sd_velocity_um_per_min=(col_vel, "std")
        )
)

# Sort by from->to
interval_stats = interval_stats.sort_values([col_from, col_to]).reset_index(drop=True)
interval_stats["interval_label"] = interval_stats[col_from].astype(str) + "→" + interval_stats[col_to].astype(str)

# Save per-interval means
interval_stats.to_csv(out_interval_csv, index=False)
print(f"Saved per-interval means CSV: {out_interval_csv}")

# Overall mean across ALL steps (not per-interval means)
overall_mean = float(work[col_vel].mean())
print(f"Overall mean velocity across ALL steps: {overall_mean:.6g} µm/min")

# IQR outlier reporting on per-step velocities (no removal)
q1 = work[col_vel].quantile(0.25)
q3 = work[col_vel].quantile(0.75)
iqr = q3 - q1
lower = q1 - 1.5 * iqr
upper = q3 + 1.5 * iqr
n_outliers = int(((work[col_vel] < lower) | (work[col_vel] > upper)).sum())
print(
    "IQR outlier report (per-step velocities; not removed): "
    f"Q1={q1:.4g}, Q3={q3:.4g}, IQR={iqr:.4g}, lower={lower:.4g}, upper={upper:.4g}, "
    f"outliers={n_outliers}/{len(work)}"
)

# -------------------- Plot --------------------
# Publication defaults
sns.set_style("whitegrid")
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

# Use numeric x (from frame) with jitter; one dot per interval
plot_df = interval_stats.copy()
plot_df["x"] = plot_df[col_from].astype(int)

plt.figure(figsize=(10, 6))
ax = sns.stripplot(
    data=plot_df,
    x="x",
    y="mean_velocity_um_per_min",
    jitter=0.25,
    size=6,
    color="#0173B2",  # colorblind-safe blue
    edgecolor="black",
    linewidth=0.5,
)

# Overall mean line
ax.axhline(overall_mean, color="black", linewidth=2.0)

ax.set_xlabel("Frame interval start (t)")
ax.set_ylabel("Mean step velocity (µm/min)")
ax.set_title(
    "Comparison of mean velocities from image to image calculated from TrackMate tracking of GT images (blue). "
    "Each dot represents the mean velocity from one time point to the next. "
    "Black line shows overall mean velocity."
)

# Improve x tick labels: show t→t+1 where possible
# Build mapping from x to label using the most common to-frame for each from-frame
label_map = (
    plot_df.groupby("x")[col_to]
           .agg(lambda s: int(pd.Series.mode(s)[0]) if len(s) else np.nan)
           .to_dict()
)
xticks = sorted(plot_df["x"].unique())
xticklabels = [f"{t}→{label_map.get(t, t+1)}" for t in xticks]
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, rotation=45, ha="right")

sns.despine(ax=ax)
plt.tight_layout()

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_svg, bbox_inches="tight")
print(f"Saved figure: {out_png}")
print(f"Saved figure: {out_svg}")
