# Stage 1 (statistics-only): Build clean per-step velocity table from TrackMate spots export.
# Inputs:  /app/data/projects/mouse_cell_tracking_trackmate/mouse_cell_spots.csv
# Outputs: /app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv
#          /app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Summary.csv
# Notes:
# - TrackMate spots export contains 3 non-data header rows (two label rows + units row).
# - Computes per-step displacement between consecutive frames within each TRACK_ID.
# - Uses pixel_size_um=0.645 and dt=5 min/frame.

import os

# ---------------------------
# Parameters (project state)
# ---------------------------
INPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/mouse_cell_spots.csv"
OUT_RESULTS = r"/app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Results.csv"
OUT_SUMMARY = r"/app/data/projects/mouse_cell_tracking_trackmate/data/Statistics_Summary.csv"

PIXEL_SIZE_UM = 0.645
DT_MIN = 5.0

# ---------------------------
# Load & clean TrackMate export
# ---------------------------
df_raw = pd.read_csv(INPUT_CSV, dtype=str)
print(f"Loaded raw spots CSV: {INPUT_CSV}")
print(f"Raw shape: {df_raw.shape}")

required_cols = ["TRACK_ID", "FRAME", "POSITION_X", "POSITION_Y"]
missing = [c for c in required_cols if c not in df_raw.columns]
if missing:
    raise ValueError(f"Missing required columns in spots CSV: {missing}. Available columns: {list(df_raw.columns)}")

# Convert to numeric; header/unit rows become NaN and will be dropped.
for c in required_cols:
    df_raw[c] = pd.to_numeric(df_raw[c], errors="coerce")

# Drop non-data rows (3 header rows + any malformed rows)
df = df_raw.dropna(subset=required_cols).copy()

# Enforce integer types for IDs/frames after dropping NaNs
# (TrackMate exports are integer-like but stored as strings)
df["TRACK_ID"] = df["TRACK_ID"].astype(int)
df["FRAME"] = df["FRAME"].astype(int)

# Positions as float
for c in ["POSITION_X", "POSITION_Y"]:
    df[c] = df[c].astype(float)

print(f"After dropping non-data rows: {df.shape}")
print("Unique tracks:", df["TRACK_ID"].nunique())

# ---------------------------
# Compute per-step velocities
# ---------------------------
# Sort within track by frame
_df = df.sort_values(["TRACK_ID", "FRAME"]).copy()

# Next-frame values within each track
_df["FRAME_NEXT"] = _df.groupby("TRACK_ID")["FRAME"].shift(-1)
_df["X_NEXT"] = _df.groupby("TRACK_ID")["POSITION_X"].shift(-1)
_df["Y_NEXT"] = _df.groupby("TRACK_ID")["POSITION_Y"].shift(-1)

# Keep only rows that have a next spot
steps = _df.dropna(subset=["FRAME_NEXT", "X_NEXT", "Y_NEXT"]).copy()
steps["FRAME_NEXT"] = steps["FRAME_NEXT"].astype(int)

# Only consecutive frames (frame-to-frame velocities)
steps["FRAME_DELTA"] = steps["FRAME_NEXT"] - steps["FRAME"]
steps = steps.loc[steps["FRAME_DELTA"] == 1].copy()

# Compute dx, dy, distance
steps["STEP_DX_PX"] = steps["X_NEXT"] - steps["POSITION_X"]
steps["STEP_DY_PX"] = steps["Y_NEXT"] - steps["POSITION_Y"]
steps["STEP_DISTANCE_PX"] = np.sqrt(steps["STEP_DX_PX"]**2 + steps["STEP_DY_PX"]**2)

# Convert to um/min
steps["STEP_VELOCITY_UM_PER_MIN"] = steps["STEP_DISTANCE_PX"] * PIXEL_SIZE_UM / DT_MIN

# IQR outlier reporting (do not remove; just report)
vel = steps["STEP_VELOCITY_UM_PER_MIN"].astype(float)
q1, q3 = np.nanpercentile(vel, [25, 75])
iqr = q3 - q1
low = q1 - 1.5 * iqr
high = q3 + 1.5 * iqr
n_outliers = int(((vel < low) | (vel > high)).sum())
print(f"Step velocities computed (consecutive frames only): N={len(steps)}")
print(f"IQR outlier check (threshold=1.5): Q1={q1:.4f}, Q3={q3:.4f}, IQR={iqr:.4f}, low={low:.4f}, high={high:.4f}, outliers={n_outliers}")

# ---------------------------
# Save clean rectangular outputs
# ---------------------------
os.makedirs(os.path.dirname(OUT_RESULTS), exist_ok=True)

out_cols = [
    "TRACK_ID",
    "FRAME",
    "FRAME_NEXT",
    "STEP_DISTANCE_PX",
    "STEP_VELOCITY_UM_PER_MIN",
]

results = steps[out_cols].copy()
# Ensure numeric dtypes in output
results["STEP_DISTANCE_PX"] = results["STEP_DISTANCE_PX"].astype(float)
results["STEP_VELOCITY_UM_PER_MIN"] = results["STEP_VELOCITY_UM_PER_MIN"].astype(float)

results.to_csv(OUT_RESULTS, index=False)
print(f"Wrote clean per-step table: {OUT_RESULTS}")

# Summary (separate file; do not mix into results)
summary = pd.DataFrame([
    {
        "N_STEPS": int(len(results)),
        "N_TRACKS": int(df["TRACK_ID"].nunique()),
        "PIXEL_SIZE_UM": float(PIXEL_SIZE_UM),
        "DT_MIN": float(DT_MIN),
        "MEAN_STEP_VELOCITY_UM_PER_MIN": float(results["STEP_VELOCITY_UM_PER_MIN"].mean()) if len(results) else np.nan,
        "MEDIAN_STEP_VELOCITY_UM_PER_MIN": float(results["STEP_VELOCITY_UM_PER_MIN"].median()) if len(results) else np.nan,
        "SD_STEP_VELOCITY_UM_PER_MIN": float(results["STEP_VELOCITY_UM_PER_MIN"].std(ddof=1)) if len(results) > 1 else np.nan,
        "IQR_OUTLIERS_N": int(n_outliers),
        "IQR_LOW": float(low),
        "IQR_HIGH": float(high),
    }
])
summary.to_csv(OUT_SUMMARY, index=False)
print(f"Wrote summary table: {OUT_SUMMARY}")

# Explicit stdout reporting for supervisor logs
if len(results):
    print(
        "Overall mean step velocity (um/min): "
        f"{summary.loc[0, 'MEAN_STEP_VELOCITY_UM_PER_MIN']:.6f} (N steps={int(summary.loc[0, 'N_STEPS'])}, N tracks={int(summary.loc[0, 'N_TRACKS'])})"
    )
else:
    print("No consecutive-frame steps found (FRAME_NEXT-FRAME != 1). Check FRAME column and tracking export.")
