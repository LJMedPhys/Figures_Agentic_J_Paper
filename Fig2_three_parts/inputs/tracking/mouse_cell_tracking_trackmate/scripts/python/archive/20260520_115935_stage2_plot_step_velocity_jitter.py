# STAGE 2 (Plotting only): Jitter plot of per-step velocity (µm/min)
# Reads: /app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv
# Writes: /app/data/projects/mouse_cell_tracking_trackmate/figures/StepVelocity_Jitter_AllSteps.png/.svg

# NOTE: Do NOT re-import pandas/numpy/matplotlib/seaborn/scipy.stats/os; they are preloaded in the environment.

# -----------------------------
# Config
# -----------------------------
input_csv = r'/app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv'
fig_dir = r'/app/data/projects/mouse_cell_tracking_trackmate/figures'
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, 'StepVelocity_Jitter_AllSteps.png')
out_svg = os.path.join(fig_dir, 'StepVelocity_Jitter_AllSteps.svg')

# -----------------------------
# Load & tidy
# -----------------------------
# The CSV was saved with a MultiIndex and a single string column; reconstruct a tidy DataFrame.
df_raw = pd.read_csv(input_csv)

# If the file is in the expected "single string column" format, split it into columns.
if df_raw.shape[1] == 1:
    col0 = df_raw.columns[0]
    s = df_raw[col0].astype(str)

    # Split on commas into expected fields.
    parts = s.str.split(',', expand=True)
    # Expected 7 fields: TRACK_ID, from_FRAME, to_FRAME, step_distance_px, step_distance_um, velocity_um_per_min, overall_mean_velocity_um_per_min
    if parts.shape[1] < 6:
        raise ValueError(
            f'Unexpected format in {input_csv}: could not split into >=6 fields. '
            f'Got {parts.shape[1]} fields. First rows: {s.head().tolist()}'
        )

    # Assign names for at least first 6; keep any extras.
    base_cols = [
        'TRACK_ID', 'from_FRAME', 'to_FRAME',
        'step_distance_px', 'step_distance_um', 'velocity_um_per_min',
        'overall_mean_velocity_um_per_min'
    ]
    use_cols = base_cols[:parts.shape[1]]
    parts.columns = use_cols

    df = parts.copy()
else:
    # Already in wide/tidy format
    df = df_raw.copy()

# Coerce numeric columns
for c in ['velocity_um_per_min', 'overall_mean_velocity_um_per_min']:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')

vel = df['velocity_um_per_min'].dropna().astype(float)

if vel.empty:
    raise ValueError('No numeric values found in velocity_um_per_min after parsing.')

mean_vel = float(vel.mean())

print(f'N steps (non-NaN velocity_um_per_min) = {vel.shape[0]}')
print(f'Overall mean velocity (µm/min) = {mean_vel:.6g}')

# -----------------------------
# Plot
# -----------------------------
sns.set_theme(style='whitegrid')

# Publication typography
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'axes.linewidth': 1.0,
})

# Use explicit blue for dots (colorblind-safe)
blue = '#0173B2'

plot_df = pd.DataFrame({
    'Group': ['All steps'] * vel.shape[0],
    'Velocity (µm/min)': vel.values
})

plt.figure(figsize=(6, 6))
ax = sns.stripplot(
    data=plot_df,
    x='Group',
    y='Velocity (µm/min)',
    jitter=0.25,
    size=3.5,
    color=blue,
    alpha=0.75,
    edgecolor=None,
)

# Mean line (black)
ax.axhline(mean_vel, color='black', linewidth=2.0, linestyle='-')

ax.set_xlabel('')
ax.set_ylabel('Velocity (µm/min)', fontsize=14)

# Keep x tick label readable
ax.set_xticklabels(['All steps'])

# Tight layout and save
plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_svg, bbox_inches='tight')
plt.close()

print(f'Saved: {out_png}')
print(f'Saved: {out_svg}')
