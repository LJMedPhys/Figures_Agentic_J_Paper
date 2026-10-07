# TrackMate tracks CSV: add velocity columns in um/min and distance columns in um
# Inputs:
#   /app/data/projects/mouse_cell_tracking_trackmate/mouse_cell_tracks.csv
# Outputs:
#   /app/data/projects/mouse_cell_tracking_trackmate/data/Tracks_with_velocity_um_per_min.csv
#   /app/data/projects/mouse_cell_tracking_trackmate/logs/Velocity_Summary.txt
# Parameters:
#   pixel_size_um=0.645
#   frame_interval_min=5

import os

INPUT_CSV = r"/app/data/projects/mouse_cell_tracking_trackmate/mouse_cell_tracks.csv"
PROJECT_DIR = r"/app/data/projects/mouse_cell_tracking_trackmate"
OUTPUT_CSV = os.path.join(PROJECT_DIR, "data", "Tracks_with_velocity_um_per_min.csv")
SUMMARY_TXT = os.path.join(PROJECT_DIR, "logs", "Velocity_Summary.txt")

pixel_size_um = 0.645
frame_interval_min = 5

# Conversion: (pixel/frame) * (um/pixel) / (min/frame) = um/min
speed_factor = pixel_size_um / frame_interval_min
# Conversion: pixel * (um/pixel) = um
length_factor = pixel_size_um

# Ensure output dirs exist
os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
os.makedirs(os.path.dirname(SUMMARY_TXT), exist_ok=True)

# Read as strings first (TrackMate exports often include header/unit rows)
df_raw = pd.read_csv(INPUT_CSV, dtype=str)

# Coerce numeric columns; non-numeric rows become NaN and will be dropped
numeric_cols = [
    "TRACK_MEAN_SPEED", "TRACK_MAX_SPEED", "TRACK_MIN_SPEED", "TRACK_MEDIAN_SPEED", "TRACK_STD_SPEED",
    "MEAN_STRAIGHT_LINE_SPEED",
    "TRACK_DISPLACEMENT", "TOTAL_DISTANCE_TRAVELED", "MAX_DISTANCE_TRAVELED",
    "TRACK_DURATION", "TRACK_START", "TRACK_STOP",
]

missing = [c for c in numeric_cols if c not in df_raw.columns]
if missing:
    raise ValueError(f"Missing expected columns in input CSV: {missing}. Found columns: {list(df_raw.columns)}")

# Convert numeric columns
for c in numeric_cols:
    df_raw[c] = pd.to_numeric(df_raw[c], errors="coerce")

# Identify and drop non-data rows (e.g., header/units rows)
# Use TRACK_MEAN_SPEED as primary indicator; require at least one key numeric field present
before_n = len(df_raw)
df = df_raw.dropna(subset=["TRACK_MEAN_SPEED", "TRACK_DISPLACEMENT"], how="any").copy()
after_n = len(df)
print(f"Rows read: {before_n}; data rows retained after dropping non-numeric header/unit rows: {after_n}")

# Preserve identifiers; ensure they remain as strings
id_cols = [c for c in ["TRACK_ID", "TRACK_INDEX", "LABEL"] if c in df.columns]
for c in id_cols:
    df[c] = df[c].astype(str)

# Add converted speed columns (um/min)
df["TRACK_MEAN_SPEED_um_per_min"] = df["TRACK_MEAN_SPEED"] * speed_factor
df["TRACK_MAX_SPEED_um_per_min"] = df["TRACK_MAX_SPEED"] * speed_factor
df["TRACK_MIN_SPEED_um_per_min"] = df["TRACK_MIN_SPEED"] * speed_factor
df["TRACK_MEDIAN_SPEED_um_per_min"] = df["TRACK_MEDIAN_SPEED"] * speed_factor
df["TRACK_STD_SPEED_um_per_min"] = df["TRACK_STD_SPEED"] * speed_factor
df["MEAN_STRAIGHT_LINE_SPEED_um_per_min"] = df["MEAN_STRAIGHT_LINE_SPEED"] * speed_factor

# Convert distance-like columns (um)
df["TRACK_DISPLACEMENT_um"] = df["TRACK_DISPLACEMENT"] * length_factor
df["TOTAL_DISTANCE_TRAVELED_um"] = df["TOTAL_DISTANCE_TRAVELED"] * length_factor
df["MAX_DISTANCE_TRAVELED_um"] = df["MAX_DISTANCE_TRAVELED"] * length_factor

# Save augmented CSV
# Keep original columns plus new ones; do not reorder aggressively
new_cols = [
    "TRACK_MEAN_SPEED_um_per_min", "TRACK_MAX_SPEED_um_per_min", "TRACK_MIN_SPEED_um_per_min",
    "TRACK_MEDIAN_SPEED_um_per_min", "TRACK_STD_SPEED_um_per_min", "MEAN_STRAIGHT_LINE_SPEED_um_per_min",
    "TRACK_DISPLACEMENT_um", "TOTAL_DISTANCE_TRAVELED_um", "MAX_DISTANCE_TRAVELED_um",
]

for c in new_cols:
    if c not in df.columns:
        raise RuntimeError(f"Failed to create expected output column: {c}")

# Write CSV
# Use float formatting for readability
float_cols = [c for c in df.columns if c.endswith("_um") or c.endswith("_um_per_min")]
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved: {OUTPUT_CSV}")

# Summary stats on mean speeds (um/min)
mean_speeds = df["TRACK_MEAN_SPEED_um_per_min"].dropna().astype(float)
summary = {
    "n_tracks": int(mean_speeds.shape[0]),
    "min_mean_speed_um_per_min": float(mean_speeds.min()) if len(mean_speeds) else float("nan"),
    "median_mean_speed_um_per_min": float(mean_speeds.median()) if len(mean_speeds) else float("nan"),
    "mean_mean_speed_um_per_min": float(mean_speeds.mean()) if len(mean_speeds) else float("nan"),
}

summary_lines = [
    "Velocity summary (TrackMate tracks)",
    f"pixel_size_um = {pixel_size_um}",
    f"frame_interval_min = {frame_interval_min}",
    f"speed_factor (um/min per pixel/frame) = {speed_factor}",
    "",
    f"N tracks (with numeric mean speed) = {summary['n_tracks']}",
    f"Min of track mean speeds (um/min) = {summary['min_mean_speed_um_per_min']:.6g}",
    f"Median of track mean speeds (um/min) = {summary['median_mean_speed_um_per_min']:.6g}",
    f"Mean of track mean speeds (um/min) = {summary['mean_mean_speed_um_per_min']:.6g}",
]

with open(SUMMARY_TXT, "w", encoding="utf-8") as f:
    f.write("\n".join(summary_lines) + "\n")

print(f"Saved: {SUMMARY_TXT}")
print("Summary:")
print("\n".join(summary_lines))
