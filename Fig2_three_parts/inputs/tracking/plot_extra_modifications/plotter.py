import os 
import numpy as np 
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 

# ----------------------------
# Paths
# ----------------------------
PROJECT_DIR = "plot_extra_modification/"
INPUT_CSV = "interval_mean_velocities.csv"
FIG_DIR = os.path.join(PROJECT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)
 
OUT_PNG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_padded_ylabel28_y0_tick1_meansegment.png")
OUT_SVG = os.path.join(FIG_DIR, "interval_mean_velocity_beeswarm_padded_ylabel28_y0_tick1_meansegment.svg")

# ----------------------------
# Load data
# ----------------------------
df = pd.read_csv(INPUT_CSV)

required_cols = ["mean_velocity_um_per_min"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}. Found columns: {list(df.columns)}")

vel = pd.to_numeric(df["mean_velocity_um_per_min"], errors="coerce")
vel = vel.replace([np.inf, -np.inf], np.nan).dropna()

if vel.empty:
    raise ValueError("No valid numeric values found in mean_velocity_um_per_min after cleaning.")

overall_mean = float(vel.mean())

# ----------------------------
# Outlier reporting (IQR) - report only, do not remove
# ----------------------------
q1 = float(vel.quantile(0.25))
q3 = float(vel.quantile(0.75))
iqr = q3 - q1
if iqr == 0:
    outlier_mask = pd.Series([False] * len(vel), index=vel.index)
else:
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    outlier_mask = (vel < lower) | (vel > upper)

n_outliers = int(outlier_mask.sum())
print(f"N intervals (points) = {len(vel)}")
print(f"Overall mean(mean_velocity_um_per_min) = {overall_mean:.6f} µm/min")
print(f"IQR outliers (threshold=1.5): {n_outliers} / {len(vel)}")

# ----------------------------
# Plot styling
# ----------------------------
sns.set_style("white")
sns.set_palette("colorblind")

TARGET_W_PX, TARGET_H_PX = 530, 597
DPI = 300
figsize = (TARGET_W_PX / DPI, TARGET_H_PX / DPI)  # (1.767, 1.990) inches

plt.rcParams.update({
    "font.size": 7,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    # Tick marks
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

fig, ax = plt.subplots(figsize=figsize, dpi=DPI)

plot_df = pd.DataFrame({
    "group": [""] * len(vel),
    "mean_velocity_um_per_min": vel.values
})

sns.swarmplot(
    data=plot_df,
    x="group",
    y="mean_velocity_um_per_min",
    ax=ax,
    size=3,           # smaller dots to fit the narrower figure
    color="#0173B2",
    edgecolor=None,
    linewidth=0,
)

# Mean segment
segment_half_width = 0.18
ax.hlines(
    overall_mean,
    xmin=0 - segment_half_width,
    xmax=0 + segment_half_width,
    colors="black",
    linewidth=1.5,
    zorder=5,
)

# ----------------------------
# Axes formatting
# ----------------------------
ax.set_xlabel("")
ax.set_ylabel("Mean velocity (µm/min)", fontsize=9)
ax.set_title("")

ax.set_ylim(bottom=0)
upper_tick = 3
yticks = np.arange(0, upper_tick + 1.0, 1.0)
ax.set_yticks(yticks)

# Force y tick marks explicitly — seaborn can reset these
ax.tick_params(
    axis="y",
    which="major",
    length=4,
    width=0.8,
    direction="out",
    left=True,       # explicitly show tick marks on left spine
    labelleft=True,
)

# Remove x ticks/labels entirely
ax.tick_params(
    axis="x",
    which="both",
    bottom=False,
    labelbottom=False,
)
ax.set_xticks([])

# Minimal spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.8)
ax.spines["bottom"].set_linewidth(0.8)

ax.grid(False)
ax.set_xlim(-0.6, 0.6)

fig.subplots_adjust(left=0.38, right=0.92, bottom=0.06, top=0.98)

# Save WITHOUT bbox_inches="tight" so tick marks aren't clipped
fig.savefig(OUT_PNG, dpi=DPI)
fig.savefig(OUT_SVG)
plt.close(fig)
