# Stage 1 (statistics/prep): compute frame-to-frame step velocities from TrackMate Spots CSV
# Environment note: pd, np, stats, os are pre-imported in the execution environment.

# -----------------------------
# User parameters (from PROJECT STATE)
# -----------------------------
INPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/mouse_cell_spots.csv"
OUTPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv"

PIXEL_SIZE_UM = 0.645
FRAME_INTERVAL_MIN = 5.0

# -----------------------------
# Helpers
# -----------------------------
def _to_numeric(series):
    """Coerce to numeric safely; non-parsable values become NaN."""
    return pd.to_numeric(series, errors='coerce')


def _read_trackmate_spots_skip_unit_rows(path):
    """Read TrackMate spots CSV and robustly skip the first header/unit rows.

    TrackMate exports often contain 3 non-data rows:
      row0: verbose header names
      row1: short header names
      row2: units in parentheses
    followed by numeric data.

    We read with header=0 to get canonical column names, then drop rows where
    TRACK_ID or POSITION_X/Y are non-numeric.
    """
    df = pd.read_csv(path, dtype=str)

    # Standardize column names (strip whitespace)
    df.columns = [c.strip() for c in df.columns]

    required = ["TRACK_ID", "POSITION_X", "POSITION_Y"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in spots CSV: {missing}. Found columns: {list(df.columns)}")

    # Coerce key columns to numeric to identify real data rows
    track_id_num = _to_numeric(df["TRACK_ID"])
    x_num = _to_numeric(df["POSITION_X"])
    y_num = _to_numeric(df["POSITION_Y"])

    mask_data = track_id_num.notna() & x_num.notna() & y_num.notna()
    df_data = df.loc[mask_data].copy()

    # Convert key columns to numeric
    df_data["TRACK_ID"] = track_id_num.loc[mask_data].astype(int)
    df_data["POSITION_X"] = x_num.loc[mask_data].astype(float)
    df_data["POSITION_Y"] = y_num.loc[mask_data].astype(float)

    # Time columns: prefer FRAME if present and numeric; else use POSITION_T
    if "FRAME" in df_data.columns:
        df_data["FRAME"] = _to_numeric(df_data["FRAME"])
    else:
        df_data["FRAME"] = np.nan

    if "POSITION_T" in df_data.columns:
        df_data["POSITION_T"] = _to_numeric(df_data["POSITION_T"])
    else:
        df_data["POSITION_T"] = np.nan

    # If FRAME is missing/NaN, fall back to POSITION_T
    df_data["FRAME_OR_T"] = df_data["FRAME"]
    df_data.loc[df_data["FRAME_OR_T"].isna(), "FRAME_OR_T"] = df_data.loc[df_data["FRAME_OR_T"].isna(), "POSITION_T"]

    # Drop rows still missing time
    before = len(df_data)
    df_data = df_data.loc[df_data["FRAME_OR_T"].notna()].copy()
    after = len(df_data)
    if after < before:
        print(f"Dropped {before-after} rows with non-numeric FRAME and POSITION_T.")

    # Ensure integer-like frames where possible
    # (TrackMate stores frames as integers; POSITION_T may be float but in frames here)
    df_data["FRAME_OR_T"] = df_data["FRAME_OR_T"].astype(float)

    return df_data


def _iqr_outlier_count(x, k=1.5):
    x = pd.Series(x).dropna().astype(float)
    if len(x) < 4:
        return 0, np.nan, np.nan, np.nan, np.nan
    q1 = x.quantile(0.25)
    q3 = x.quantile(0.75)
    iqr = q3 - q1
    lo = q1 - k * iqr
    hi = q3 + k * iqr
    n_out = int(((x < lo) | (x > hi)).sum())
    return n_out, float(q1), float(q3), float(lo), float(hi)


# -----------------------------
# Main
# -----------------------------
spots = _read_trackmate_spots_skip_unit_rows(INPUT_CSV)

# Sort within track by time
spots = spots.sort_values(["TRACK_ID", "FRAME_OR_T", "ID" if "ID" in spots.columns else "POSITION_X"]).reset_index(drop=True)

# Compute step-wise differences within each track
spots["from_FRAME"] = spots.groupby("TRACK_ID")["FRAME_OR_T"].shift(0)
spots["to_FRAME"] = spots.groupby("TRACK_ID")["FRAME_OR_T"].shift(-1)
spots["x0"] = spots.groupby("TRACK_ID")["POSITION_X"].shift(0)
spots["y0"] = spots.groupby("TRACK_ID")["POSITION_Y"].shift(0)
spots["x1"] = spots.groupby("TRACK_ID")["POSITION_X"].shift(-1)
spots["y1"] = spots.groupby("TRACK_ID")["POSITION_Y"].shift(-1)

steps = spots.loc[spots["to_FRAME"].notna()].copy()

# Only keep consecutive detections (t -> t+1) in frame units
steps["delta_frame"] = steps["to_FRAME"] - steps["from_FRAME"]
steps = steps.loc[np.isclose(steps["delta_frame"].astype(float), 1.0)].copy()

# Step displacement in pixels
steps["dx_px"] = steps["x1"] - steps["x0"]
steps["dy_px"] = steps["y1"] - steps["y0"]
steps["step_distance_px"] = np.sqrt(steps["dx_px"]**2 + steps["dy_px"]**2)

# Convert to um and velocity um/min
steps["step_distance_um"] = steps["step_distance_px"] * PIXEL_SIZE_UM
steps["velocity_um_per_min"] = steps["step_distance_um"] / FRAME_INTERVAL_MIN

# Long-format output table
long_cols = [
    "TRACK_ID",
    "from_FRAME",
    "to_FRAME",
    "step_distance_px",
    "step_distance_um",
    "velocity_um_per_min",
]
long_df = steps[long_cols].copy()

# Overall mean velocity across all steps
overall_mean_vel = float(long_df["velocity_um_per_min"].mean()) if len(long_df) else np.nan
long_df["overall_mean_velocity_um_per_min"] = overall_mean_vel

# Descriptive stats
v = long_df["velocity_um_per_min"].dropna().astype(float)
summary = {
    "metric": "velocity_um_per_min",
    "N_steps": int(v.shape[0]),
    "mean": float(v.mean()) if len(v) else np.nan,
    "median": float(v.median()) if len(v) else np.nan,
    "std": float(v.std(ddof=1)) if len(v) > 1 else np.nan,
    "min": float(v.min()) if len(v) else np.nan,
    "max": float(v.max()) if len(v) else np.nan,
    "pixel_size_um": PIXEL_SIZE_UM,
    "frame_interval_min": FRAME_INTERVAL_MIN,
}
summary_df = pd.DataFrame([summary])

# Outlier reporting (IQR method) - report only, do not remove
n_out, q1, q3, lo, hi = _iqr_outlier_count(v, k=1.5)
print(f"Computed {int(v.shape[0])} step velocities from {long_df['TRACK_ID'].nunique()} tracks (consecutive frames only).")
print(f"Overall mean velocity (um/min): {overall_mean_vel}")
print(f"IQR outlier detection on velocity_um_per_min (k=1.5): n_outliers={n_out}, Q1={q1}, Q3={q3}, lo={lo}, hi={hi}")

outlier_df = pd.DataFrame([{
    "metric": "velocity_um_per_min",
    "iqr_k": 1.5,
    "n_outliers": n_out,
    "Q1": q1,
    "Q3": q3,
    "lower_bound": lo,
    "upper_bound": hi,
}])

# Write a single CSV with two blocks separated by a blank line and section headers.
# This is robust for human reading and still parseable by downstream scripts.
os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
with open(OUTPUT_CSV, "w", newline="") as f:
    f.write("# LongFormat_StepVelocities\n")
    long_df.to_csv(f, index=False)
    f.write("\n# Summary_DescriptiveStats\n")
    summary_df.to_csv(f, index=False)
    f.write("\n# Outliers_IQR_Report\n")
    outlier_df.to_csv(f, index=False)

print(f"Saved Statistics_Results.csv to: {OUTPUT_CSV}")
