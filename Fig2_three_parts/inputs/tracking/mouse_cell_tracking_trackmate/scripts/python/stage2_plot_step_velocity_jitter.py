# Stage 2 (plotting-only): Publication-quality jitter/strip plot of per-step velocities.
# Reads per-step STEP_VELOCITY_UM_PER_MIN from Statistics_Results.csv, computes overall mean,
# plots jittered points (single category) with mean line, and saves mean to CSV.
# Environment presets available: pd, np, plt, sns, stats, os

input_csv = r'/app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv'
project_dir = r'/app/data/projects/mouse_cell_tracking_trackmate'
fig_dir = os.path.join(project_dir, 'figures')
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, 'step_velocity_jitter.png')
out_svg = os.path.join(fig_dir, 'step_velocity_jitter.svg')
out_mean_csv = os.path.join(fig_dir, 'step_velocity_mean_used.csv')


def _read_step_velocity_column_robust(path):
    """Robustly read STEP_VELOCITY_UM_PER_MIN from a potentially messy/mixed CSV.

    Primary path: normal pd.read_csv.
    Fallback: python engine with on_bad_lines='skip'.

    Returns
    -------
    vel : pd.Series (float) of non-NaN velocities.
    """
    # Attempt 1: standard read
    try:
        df = pd.read_csv(path)
    except Exception:
        # Attempt 2: tolerate malformed lines
        try:
            df = pd.read_csv(path, engine='python', on_bad_lines='skip')
        except TypeError:
            # Older pandas compatibility
            df = pd.read_csv(path, engine='python', error_bad_lines=False, warn_bad_lines=True)

    df.columns = [str(c).strip() for c in df.columns]

    target = 'STEP_VELOCITY_UM_PER_MIN'
    if target not in df.columns:
        # Try case-insensitive match
        lower_map = {str(c).strip().lower(): c for c in df.columns}
        if target.lower() in lower_map:
            target = lower_map[target.lower()]
        else:
            raise ValueError(
                f"Could not find STEP_VELOCITY_UM_PER_MIN column. Columns seen: {list(df.columns)}"
            )

    vel = pd.to_numeric(df[target], errors='coerce')
    vel = vel.loc[vel.notna()].astype(float)
    return vel


vel = _read_step_velocity_column_robust(input_csv)

n_steps = int(vel.shape[0])
mean_vel = float(vel.mean())

print(f"Loaded step velocities: n={n_steps}")
print(f"Mean STEP_VELOCITY_UM_PER_MIN used for mean line: {mean_vel:.6g} µm/min")

# Save mean used for the line
mean_df = pd.DataFrame([
    {
        'metric': 'STEP_VELOCITY_UM_PER_MIN',
        'n_steps': n_steps,
        'mean_um_per_min': mean_vel,
        'source_csv': os.path.basename(input_csv)
    }
])
mean_df.to_csv(out_mean_csv, index=False)
print(f"Saved mean CSV: {out_mean_csv}")

# Publication styling
sns.set_theme(style='whitegrid')
sns.set_palette('colorblind')
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'axes.linewidth': 1.0,
    'lines.linewidth': 1.5,
})

plot_df = pd.DataFrame({
    'Category': ['Velocity'] * n_steps,
    'STEP_VELOCITY_UM_PER_MIN': vel.values
})

plt.figure(figsize=(8, 6))
ax = sns.stripplot(
    data=plot_df,
    x='Category',
    y='STEP_VELOCITY_UM_PER_MIN',
    jitter=0.25,
    size=4,
    alpha=0.8,
    color='#0173B2',  # colorblind-safe blue
    edgecolor='none'
)

# Mean line
ax.axhline(mean_vel, color='black', linewidth=2.0, zorder=3)

ax.set_xlabel('')
ax.set_ylabel('Step velocity (µm/min)', fontsize=14)
ax.set_title('TrackMate step velocities (all consecutive-frame steps)', fontsize=16)

sns.despine(ax=ax)

plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_svg, bbox_inches='tight')
plt.close()

print(f"Saved figure: {out_png}")
print(f"Saved figure: {out_svg}")
