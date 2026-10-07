# Stage 2: Publication-quality plotting (NO statistics computed here)
# Robustly loads only long-format per-step velocities from a mixed-content CSV.

# Environment presets available: pd, np, plt, sns, stats, os

input_csv = r'/app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv'
project_dir = r'/app/data/projects/mouse_cell_tracking_trackmate'
fig_dir = os.path.join(project_dir, 'figures')
os.makedirs(fig_dir, exist_ok=True)

out_png = os.path.join(fig_dir, 'step_velocity_jitter.png')
out_svg = os.path.join(fig_dir, 'step_velocity_jitter.svg')


def _read_mixed_csv_longformat(path):
    """Read a CSV that contains a long-format block plus a summary block with different columns.

    Strategy:
      1) Read with python engine and on_bad_lines='skip' to bypass malformed rows.
      2) If the file is actually a single-column dump (common when commas/quotes are inconsistent),
         split that column by comma.
      3) Coerce velocity_um_per_min to numeric and keep only numeric rows.

    Returns
    -------
    df_long : pd.DataFrame with at least ['velocity_um_per_min'] numeric.
    """
    # Attempt 1: normal read, skipping bad lines
    try:
        df = pd.read_csv(path, engine='python', on_bad_lines='skip')
    except TypeError:
        # Older pandas compatibility
        df = pd.read_csv(path, engine='python', error_bad_lines=False, warn_bad_lines=True)

    # If the CSV got parsed into a single column, try to split it
    if df.shape[1] == 1:
        only_col = df.columns[0]
        # Split into columns by comma; keep as strings
        split_df = df[only_col].astype(str).str.split(',', expand=True)
        # Use first row as header if it looks like it contains the velocity column name
        header_row = split_df.iloc[0].astype(str).tolist()
        if any('velocity_um_per_min' in str(x) for x in header_row):
            split_df.columns = header_row
            split_df = split_df.iloc[1:].copy()
            df = split_df
        else:
            # If no header detected, keep original df; downstream will fail with clear message
            df = split_df

    # Normalize column names (strip whitespace)
    df.columns = [str(c).strip() for c in df.columns]

    if 'velocity_um_per_min' not in df.columns:
        # Sometimes header may have been shifted; try to find a column containing the name
        candidates = [c for c in df.columns if 'velocity' in c and 'um' in c]
        raise ValueError(
            f"Could not find 'velocity_um_per_min' column. Columns seen: {list(df.columns)[:30]}" +
            (f"; candidates={candidates}" if candidates else '')
        )

    df['velocity_um_per_min'] = pd.to_numeric(df['velocity_um_per_min'], errors='coerce')
    df_long = df.loc[df['velocity_um_per_min'].notna(), ['velocity_um_per_min']].copy()

    return df_long


df_long = _read_mixed_csv_longformat(input_csv)

n_steps = int(df_long.shape[0])
mean_vel = float(df_long['velocity_um_per_min'].mean())

print(f"Loaded per-step velocities: n={n_steps}")
print(f"Mean velocity (um/min) of plotted steps: {mean_vel:.6g}")

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

# Prepare plotting dataframe with a single category
plot_df = df_long.copy()
plot_df['Group'] = 'All steps'

plt.figure(figsize=(8, 6))
ax = sns.stripplot(
    data=plot_df,
    x='Group',
    y='velocity_um_per_min',
    jitter=0.25,
    size=4,
    alpha=0.8,
    color='#0173B2',  # colorblind-safe blue
    edgecolor='none'
)

# Mean line (horizontal) across the category
ax.axhline(mean_vel, color='black', linewidth=2.0, zorder=3)

ax.set_xlabel('')
ax.set_ylabel('Step velocity (µm/min)', fontsize=14)
ax.set_title('TrackMate step velocities (all consecutive-frame steps)', fontsize=16)

# Clean spines
sns.despine(ax=ax)

plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_svg, bbox_inches='tight')
plt.close()

print(f"Saved: {out_png}")
print(f"Saved: {out_svg}")
