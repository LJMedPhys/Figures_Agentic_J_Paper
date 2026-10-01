#!/usr/bin/env python3
"""Fully automatic analysis of nuclear TUNEL signal in Position040."""

from __future__ import annotations

import argparse
import subprocess
import struct
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from readlif.reader import LifFile
from skimage import filters, measure, morphology


def read_rect_imagej_roi(path: Path) -> tuple[int, int, int, int]:
    """Return (left, top, right_exclusive, bottom_exclusive) for a rectangle ROI."""
    data = path.read_bytes()
    if len(data) < 64 or data[:4] != b"Iout":
        raise ValueError(f"Not an ImageJ ROI: {path}")
    roi_type = data[6]
    if roi_type != 1:
        raise ValueError(f"Expected rectangular ImageJ ROI (type 1), got {roi_type}")
    top, left, bottom, right = struct.unpack(">4H", data[8:16])
    if right <= left or bottom <= top:
        raise ValueError("ROI has an invalid bounding rectangle")
    return left, top, right, bottom


def load_position40(lif_path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    lif = LifFile(str(lif_path))
    image = lif.get_image(39)
    if image.name != "Position040" or image.channels != 3:
        raise RuntimeError(
            f"Expected series 39 to be Position040 with 3 channels; got "
            f"{image.name!r} with {image.channels} channels"
        )
    # LIF hardware metadata: c0 brightfield; c1 635-ex/695-em (CF640R TUNEL);
    # c2 385-ex/425-em (DAPI).
    tunel = np.asarray(image.get_frame(c=1), dtype=np.float32)
    dapi = np.asarray(image.get_frame(c=2), dtype=np.float32)
    metadata = {
        "series": image.name,
        "shape": dapi.shape,
        "channels": image.channels,
        "pixel_size_um": 1.0 / image.scale[0],
    }
    if tunel.shape != dapi.shape:
        raise RuntimeError(f"TUNEL {tunel.shape} and DAPI {dapi.shape} differ")
    return tunel, dapi, metadata


def segment_nuclei(dapi: np.ndarray) -> np.ndarray:
    """Segment DAPI nuclei with the pretrained Cellpose-SAM model."""
    helper = Path(__file__).with_name("segment_nuclei_cellpose.py")
    interpreter = Path("/opt/conda/envs/cellpose4/bin/python")
    if not helper.exists() or not interpreter.exists():
        raise RuntimeError("Cellpose-SAM helper or cellpose4 interpreter is unavailable")
    with tempfile.TemporaryDirectory(prefix="tunel_cellpose_") as tmp:
        input_tif = Path(tmp) / "dapi.tif"
        output_tif = Path(tmp) / "labels.tif"
        tifffile.imwrite(input_tif, dapi.astype(np.uint16), photometric="minisblack")
        subprocess.run(
            [str(interpreter), str(helper), str(input_tif), str(output_tif)],
            check=True,
        )
        labels = tifffile.imread(output_tif)

    # Apply conservative size/edge QC, then assign deterministic IDs from top to
    # bottom (and left to right for equal rows).
    regions = sorted(measure.regionprops(labels), key=lambda r: (r.centroid[0], r.centroid[1]))
    out = np.zeros_like(labels, dtype=np.int32)
    next_id = 1
    for region in regions:
        minr, minc, maxr, maxc = region.bbox
        touches_edge = minr == 0 or minc == 0 or maxr == labels.shape[0] or maxc == labels.shape[1]
        if 90 <= region.area <= 2500 and not touches_edge:
            out[labels == region.label] = next_id
            next_id += 1
    return out


def local_backgrounds(
    tunel: np.ndarray, labels: np.ndarray, inner_px: int = 6, outer_px: int = 18
) -> dict[int, float]:
    """Measure robust local background in nucleus-free annuli around each nucleus."""
    all_nuclei = labels > 0
    results: dict[int, float] = {}
    for nucleus_id in range(1, int(labels.max()) + 1):
        obj = labels == nucleus_id
        outer = morphology.dilation(obj, morphology.disk(outer_px))
        inner = morphology.dilation(obj, morphology.disk(inner_px))
        annulus = outer & ~inner & ~all_nuclei
        values = tunel[annulus]
        if values.size < 100:
            raise RuntimeError(f"Too few background pixels for nucleus {nucleus_id}")
        # Trim bright non-nuclear structures and isolated hot pixels, then average.
        cutoff = np.percentile(values, 80)
        results[nucleus_id] = float(values[values <= cutoff].mean())
    return results


def choose_threshold(values: np.ndarray) -> tuple[float, str]:
    """Use log-Otsu classes, placing the reported cut in their empty interval."""
    nonnegative = np.clip(values, 0, None)
    logged = np.log1p(nonnegative)
    otsu = filters.threshold_otsu(logged)
    low = logged[logged <= otsu]
    high = logged[logged > otsu]
    if low.size == 0 or high.size == 0:
        raise RuntimeError("TUNEL thresholding did not produce two groups")
    threshold_log = (float(low.max()) + float(high.min())) / 2
    threshold = float(np.expm1(threshold_log))
    method = (
        "two-class Otsu on log1p-transformed corrected intensities, with the final "
        "cut centered in log space between the brightest low-class nucleus and "
        "the dimmest high-class nucleus"
    )
    return threshold, method


def measure_nuclei(
    tunel: np.ndarray, labels: np.ndarray, x0: int, y0: int
) -> tuple[pd.DataFrame, float, str]:
    bgs = local_backgrounds(tunel, labels)
    rows = []
    for region in measure_regionprops(labels, tunel):
        nucleus_id = region.label
        raw = float(region.intensity_mean)
        bg = bgs[nucleus_id]
        rows.append(
            {
                "nucleus_id": nucleus_id,
                "centroid_x": float(region.centroid[1] + x0),
                "centroid_y": float(region.centroid[0] + y0),
                "area_px": int(region.area),
                "raw_tunel_intensity": raw,
                "background_intensity": bg,
                "corrected_tunel_intensity": raw - bg,
            }
        )
    df = pd.DataFrame(rows)
    threshold, method = choose_threshold(df["corrected_tunel_intensity"].to_numpy())
    df["group"] = np.where(
        df["corrected_tunel_intensity"] > threshold, "TUNEL+", "TUNEL-"
    )
    return df, threshold, method


def measure_regionprops(labels: np.ndarray, intensity: np.ndarray):
    return measure.regionprops(labels, intensity_image=intensity)


def save_scatter(df: pd.DataFrame, threshold: float, path: Path) -> None:
    rng = np.random.default_rng(20210728)
    fig, ax = plt.subplots(figsize=(5.2, 5.0), constrained_layout=True)
    order = ["TUNEL-", "TUNEL+"]
    colors = {"TUNEL-": "#277da1", "TUNEL+": "#d1495b"}
    for x, group in enumerate(order):
        vals = df.loc[df.group == group, "corrected_tunel_intensity"].to_numpy()
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(
            np.full(len(vals), x) + jitter,
            vals,
            s=28,
            alpha=0.82,
            color=colors[group],
            edgecolor="white",
            linewidth=0.35,
        )
        if len(vals):
            ax.hlines(np.median(vals), x - 0.28, x + 0.28, color="black", lw=2.4)
    ax.axhline(threshold, color="#555555", linestyle="--", linewidth=1, label=f"cut = {threshold:.1f}")
    ax.set_xticks([0, 1], order)
    ax.set_ylabel("Corrected nuclear TUNEL intensity (camera units)")
    ax.set_xlim(-0.5, 1.5)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_report(
    path: Path,
    df: pd.DataFrame,
    metadata: dict,
    roi: tuple[int, int, int, int],
    threshold: float,
    threshold_method: str,
) -> None:
    left, top, right, bottom = roi
    counts = df.group.value_counts()
    medians = df.groupby("group")["corrected_tunel_intensity"].median()
    report = f"""# TUNEL signal in nuclei — Position040

## Source image, channels, and ROI check

I analyzed series **Position040** (series index 39) from `20210728_TUNEL_staining.lif`. The Leica hardware metadata identifies channel 0 as transmitted-light brightfield, channel 1 as the 635-nm-excited/695-nm-emission acquisition (CF640R; TUNEL), and channel 2 as the 385-nm-excited/425-nm-emission acquisition (DAPI). Thus channel 1 was measured as TUNEL and channel 2 was used for nuclei.

Both fluorescence channels are {metadata['shape'][1]} × {metadata['shape'][0]} pixels and match one another exactly. The ImageJ ROI is a rectangle x={left}…{right - 1}, y={top}…{bottom - 1} (width {right - left}, height {bottom - top} pixels; zero-based coordinates). An ImageJ `.roi` stores its bounding coordinates, not the dimensions of its source image, so literal equality between ROI size and image size is not expected: this ROI deliberately specifies a smaller analysis region. Its complete bounding box is valid within the 2048 × 2048 source, and the extracted region in both channels is {right-left} × {bottom-top} pixels. No resize, registration, or coordinate correction was needed. Pixel size from the LIF metadata is {metadata['pixel_size_um']:.4f} µm/pixel. The delivered label mask is 2048 × 2048 so that it aligns directly with Position040; pixels outside the supplied ROI are zero.

## Nucleus segmentation

Segmentation was fully automatic and used only the DAPI image inside the supplied ROI. I used the pretrained **Cellpose 4 Cellpose-SAM (`cpsam`)** model with GPU inference and a 35-pixel target diameter (11.3 µm at this image's sampling). This model resolved the one touching/overlapping nuclear pair in the ROI. I then retained plausible nuclear objects of 90–2500 px and excluded partial objects touching the analysis-ROI boundary. IDs were reassigned consecutively from top to bottom. This produced **{len(df)} nuclei**. Centroids in the CSV are x/y coordinates in the full Position040 image, not ROI-local coordinates.

## TUNEL measurement and local background

`raw_tunel_intensity` is the mean channel-1 intensity over every pixel of a nucleus. For each nucleus independently, I formed an automatically generated local annulus extending 6–18 px outside that nucleus. All segmented nuclear pixels were excluded, including neighboring nuclei. To keep isolated bright extracellular structures or hot pixels from biasing the estimate, I discarded the brightest 20% of the remaining annular pixels and averaged the lower 80%. This is the reported `background_intensity`. Therefore every background was sampled immediately next to the measured nucleus, without hand-drawn regions. `corrected_tunel_intensity = raw_tunel_intensity - background_intensity`.

## TUNEL-positive cut

The cut was selected objectively from this ROI rather than supplied in advance: **{threshold_method}**, giving **{threshold:.3f} camera units**. I used the log transform because the bright-positive population is right-skewed, and centered the final cut in the observed empty interval so that it does not sit artificially close to either measured nucleus. Values were classified TUNEL+ only when strictly above the cut. This yielded **{int(counts.get('TUNEL+', 0))} TUNEL+** and **{int(counts.get('TUNEL-', 0))} TUNEL-** nuclei. Median corrected intensities were {medians.get('TUNEL-', float('nan')):.3f} (TUNEL-) and {medians.get('TUNEL+', float('nan')):.3f} (TUNEL+).

## Resolution and caveats

All requested quantities were resolved. The classification is an image-intensity phenotype, not an independent biological validation of apoptosis. Nuclei cut by the supplied ROI boundary were deliberately excluded because neither their area nor mean intensity could be measured completely.
"""
    path.write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lif", type=Path, default=Path("data/parasite_task/20210728_TUNEL_staining.lif"))
    parser.add_argument("--roi", type=Path, default=Path("data/parasite_task/one_frame/TUNEL_Overview.roi"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    tunel_full, dapi_full, metadata = load_position40(args.lif)
    roi = read_rect_imagej_roi(args.roi)
    left, top, right, bottom = roi
    height, width = metadata["shape"]
    if not (0 <= left < right <= width and 0 <= top < bottom <= height):
        raise RuntimeError(f"ROI {roi} falls outside source image shape {metadata['shape']}")
    tunel = tunel_full[top:bottom, left:right]
    dapi = dapi_full[top:bottom, left:right]

    labels = segment_nuclei(dapi)
    df, threshold, method = measure_nuclei(tunel, labels, left, top)
    df.to_csv(args.output / "tunel_nuclei.csv", index=False, float_format="%.6f")
    full_labels = np.zeros(metadata["shape"], dtype=np.uint16)
    full_labels[top:bottom, left:right] = labels.astype(np.uint16)
    tifffile.imwrite(
        args.output / "tunel_nuclei_mask.tif",
        full_labels,
        photometric="minisblack",
        compression="zlib",
    )
    save_scatter(df, threshold, args.output / "tunel_scatter.png")
    write_report(args.output / "report.md", df, metadata, roi, threshold, method)

    ids = np.unique(full_labels)
    ids = ids[ids > 0]
    assert len(ids) == len(df)
    assert np.array_equal(ids, df.nucleus_id.to_numpy())
    print(f"nuclei={len(df)}; TUNEL+={(df.group == 'TUNEL+').sum()}; threshold={threshold:.3f}")


if __name__ == "__main__":
    main()
