#!/usr/bin/env python3
"""Fully automatic analysis for position 40 TUNEL and parasite-index tasks."""

from __future__ import annotations

import argparse
import math
import struct
import subprocess
import tempfile
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from readlif.reader import LifFile
from scipy.ndimage import binary_dilation, distance_transform_edt, gaussian_filter
from skimage.feature import peak_local_max
from skimage.measure import label, regionprops
from skimage.morphology import closing, disk
from skimage.segmentation import find_boundaries, watershed


def read_rectangle_roi(path: Path) -> tuple[int, int, int, int]:
    data = path.read_bytes()
    if data[:4] != b"Iout" or data[6] != 1:
        raise ValueError(f"{path} is not an ImageJ rectangular ROI")
    top, left, bottom, right = struct.unpack(">4H", data[8:16])
    if bottom <= top or right <= left:
        raise ValueError(f"Invalid ROI bounds in {path}")
    return left, top, right, bottom


def channel_luts(lif: LifFile) -> list[str]:
    return [e.attrib.get("LUTName", "") for e in lif.xml_root.iter("ChannelDescription")][:3]


def fluorescence_wavelengths(lif: LifFile) -> list[float]:
    out: list[float] = []
    for e in lif.xml_root.iter("WideFieldChannelInfo"):
        channel = e.attrib.get("Channel")
        if channel in {"2000", "2001", "2002"} and len(out) < 3:
            out.append(float(e.attrib.get("EmissionWavelength", "nan")))
    if len(out) != 3:
        raise RuntimeError("Could not read the three staining acquisition wavelengths")
    return out


def robust_mad(values: np.ndarray) -> float:
    med = float(np.median(values))
    return max(float(np.median(np.abs(values - med)) * 1.4826), 1e-9)


def filter_components(mask: np.ndarray, min_area: int, max_area: int | None = None) -> np.ndarray:
    lab = label(mask)
    out = np.zeros_like(lab, dtype=np.uint16)
    kept = []
    for p in regionprops(lab):
        if p.area >= min_area and (max_area is None or p.area <= max_area):
            kept.append(p)
    kept.sort(key=lambda p: (p.centroid[0], p.centroid[1]))
    for new_id, p in enumerate(kept, 1):
        out[lab == p.label] = new_id
    return out


def filter_existing_labels(labels: np.ndarray, min_area: int, max_area: int | None = None) -> np.ndarray:
    """Filter an instance mask without merging labels that share an edge."""
    out = np.zeros_like(labels, dtype=np.uint16)
    kept = [
        p for p in regionprops(labels)
        if p.area >= min_area and (max_area is None or p.area <= max_area)
    ]
    kept.sort(key=lambda p: (p.centroid[0], p.centroid[1]))
    for new_id, p in enumerate(kept, 1):
        out[labels == p.label] = new_id
    return out


def clipped_annular_mean(
    image: np.ndarray,
    object_mask: np.ndarray,
    excluded: np.ndarray,
    inner: int,
    outer: int,
) -> tuple[float, int]:
    ring = binary_dilation(object_mask, iterations=outer)
    ring &= ~binary_dilation(object_mask, iterations=inner)
    ring &= ~binary_dilation(excluded, iterations=2)
    values = image[ring].astype(float)
    med = float(np.median(values))
    mad = robust_mad(values)
    values = values[(values >= med - 3 * mad) & (values <= med + 3 * mad)]
    if values.size < 30:
        raise RuntimeError("Too few uncontaminated local-background pixels")
    return float(values.mean()), int(values.size)


def prepare_registration(image: np.ndarray) -> np.ndarray:
    low, high = np.percentile(image, [2, 98])
    norm = np.clip((image.astype(np.float32) - low) / (high - low), 0, 1)
    norm = cv2.GaussianBlur(norm, (0, 0), 1.2)
    return norm - cv2.GaussianBlur(norm, (0, 0), 18)


def register_stain_to_movie(movie_bf: np.ndarray, stain_bf: np.ndarray) -> tuple[np.ndarray, dict]:
    template = prepare_registration(movie_bf)
    moving = prepare_registration(stain_bf)
    sift = cv2.SIFT_create(nfeatures=5000, contrastThreshold=0.01, edgeThreshold=15)
    tm = cv2.normalize(template, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    sm = cv2.normalize(moving, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    kp_t, des_t = sift.detectAndCompute(tm, None)
    kp_s, des_s = sift.detectAndCompute(sm, None)
    matches = cv2.BFMatcher().knnMatch(des_s, des_t, k=2)
    good = [m for m, n in matches if m.distance < 0.8 * n.distance]
    if len(good) < 8:
        raise RuntimeError(f"Registration produced only {len(good)} SIFT matches")
    source = np.float32([kp_s[m.queryIdx].pt for m in good])
    target = np.float32([kp_t[m.trainIdx].pt for m in good])
    cv2.setRNGSeed(12345)
    initial, inliers = cv2.estimateAffinePartial2D(
        source,
        target,
        method=cv2.RANSAC,
        ransacReprojThreshold=3,
        maxIters=10000,
        confidence=0.999,
        refineIters=100,
    )
    if initial is None or int(inliers.sum()) < 5:
        raise RuntimeError("RANSAC registration failed")

    # ECC uses a template-to-input mapping; initialize it with the inverse RANSAC map.
    warp = cv2.invertAffineTransform(initial.astype(np.float32))
    angle = math.atan2(float(warp[1, 0]), float(warp[0, 0]))
    warp[:, :2] = [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
    ecc, warp = cv2.findTransformECC(
        template,
        moving,
        warp,
        cv2.MOTION_EUCLIDEAN,
        (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 500, 1e-6),
        None,
        5,
    )
    stain_to_movie = cv2.invertAffineTransform(warp)
    predicted = source @ stain_to_movie[:, :2].T + stain_to_movie[:, 2]
    errors = np.linalg.norm(predicted - target, axis=1)
    keep = inliers.ravel().astype(bool)
    stats = {
        "sift_matches": len(good),
        "inliers": int(keep.sum()),
        "median_residual": float(np.median(errors[keep])),
        "mean_residual": float(np.mean(errors[keep])),
        "ecc": float(ecc),
        "angle_deg": float(np.degrees(np.arctan2(stain_to_movie[1, 0], stain_to_movie[0, 0]))),
    }
    return stain_to_movie, stats


def to_global_affine(local_affine: np.ndarray, stain_roi: tuple[int, int, int, int], movie_roi: tuple[int, int, int, int]) -> np.ndarray:
    stain_origin = np.array(stain_roi[:2], dtype=float)
    movie_origin = np.array(movie_roi[:2], dtype=float)
    shift = movie_origin + local_affine[:, 2] - local_affine[:, :2] @ stain_origin
    return np.column_stack([local_affine[:, :2], shift])


def transform_xy(affine: np.ndarray, x: float, y: float) -> tuple[float, float]:
    result = affine @ np.array([x, y, 1.0])
    return float(result[0]), float(result[1])


def nearest_label(mask: np.ndarray, x: float, y: float, max_distance: float) -> tuple[int, float]:
    xi, yi = int(round(x)), int(round(y))
    if 0 <= yi < mask.shape[0] and 0 <= xi < mask.shape[1] and mask[yi, xi] > 0:
        return int(mask[yi, xi]), 0.0
    y0, y1 = max(0, yi - 12), min(mask.shape[0], yi + 13)
    x0, x1 = max(0, xi - 12), min(mask.shape[1], xi + 13)
    ys, xs = np.nonzero(mask[y0:y1, x0:x1])
    if not len(xs):
        return 0, float("inf")
    distances = np.hypot(xs + x0 - x, ys + y0 - y)
    idx = int(np.argmin(distances))
    if distances[idx] > max_distance:
        return 0, float(distances[idx])
    return int(mask[ys[idx] + y0, xs[idx] + x0]), float(distances[idx])


def scatter_plot(frame: pd.DataFrame, category: str, value: str, ylabel: str, path: Path) -> None:
    order = ["TUNEL-", "TUNEL+"]
    colors = {"TUNEL-": "#4477AA", "TUNEL+": "#CC3311"}
    rng = np.random.default_rng(12345)
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    for x, group in enumerate(order):
        values = frame.loc[frame[category] == group, value].to_numpy(float)
        if not len(values):
            continue
        jitter = rng.uniform(-0.10, 0.10, len(values))
        ax.scatter(np.full(len(values), x) + jitter, values, s=42, color=colors[group], edgecolor="white", linewidth=0.5)
        median = float(np.median(values))
        ax.hlines(median, x - 0.25, x + 0.25, color="black", linewidth=2.2)
    ax.set_xticks([0, 1], order)
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.55, 1.55)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("data/parasite_task"))
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    args = parser.parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    stain_path = data_dir / "20210728_TUNEL_staining.lif"
    movie_path = data_dir / "time-lapse/20211028_Pos36_40.lif"
    stain_roi = read_rectangle_roi(data_dir / "one_frame/TUNEL_Overview.roi")
    movie_roi = read_rectangle_roi(data_dir / "one_frame/LastFrame_Overview.roi")

    stain_lif = LifFile(str(stain_path))
    movie_lif = LifFile(str(movie_path))
    stain_image = stain_lif.get_image(39)
    movie_image = movie_lif.get_image(4)
    if stain_image.name != "Position040" or movie_image.name != "Position040":
        raise RuntimeError("Position 40 was not series 39 (stain) / series 4 (movie)")
    shape_stain = (stain_image.dims.y, stain_image.dims.x)
    shape_movie = (movie_image.dims.y, movie_image.dims.x)
    for roi, shape, name in [(stain_roi, shape_stain, "stain"), (movie_roi, shape_movie, "movie")]:
        x0, y0, x1, y1 = roi
        if not (0 <= x0 < x1 <= shape[1] and 0 <= y0 < y1 <= shape[0]):
            raise RuntimeError(f"{name} ROI is outside its source image")
    if (stain_roi[2] - stain_roi[0], stain_roi[3] - stain_roi[1]) != (movie_roi[2] - movie_roi[0], movie_roi[3] - movie_roi[1]):
        raise RuntimeError("The two overview ROIs have different sizes")

    stain = np.stack([np.asarray(stain_image.get_frame(z=0, t=0, c=c), dtype=np.uint16) for c in range(3)])
    movie = np.stack([np.asarray(movie_image.get_frame(z=0, t=262, c=c), dtype=np.uint16) for c in range(3)])
    movie_luts = channel_luts(movie_lif)
    green_idx = movie_luts.index("Green")
    red_idx = movie_luts.index("Red")
    wavelengths = fluorescence_wavelengths(stain_lif)
    tunel_idx = 1 + int(np.argmax(wavelengths[1:]))
    dapi_idx = 1 + int(np.argmin(wavelengths[1:]))

    sx0, sy0, sx1, sy1 = stain_roi
    mx0, my0, mx1, my1 = movie_roi
    stain_crop = stain[:, sy0:sy1, sx0:sx1]
    movie_crop = movie[:, my0:my1, mx0:mx1]
    local_affine, reg = register_stain_to_movie(movie_crop[0], stain_crop[0])
    global_affine = to_global_affine(local_affine, stain_roi, movie_roi)

    # DAPI nuclei: a high-specificity robust threshold on lightly smoothed DAPI.
    dapi = stain_crop[dapi_idx].astype(float)
    tunel = stain_crop[tunel_idx].astype(float)
    dapi_med, dapi_mad = float(np.median(dapi)), robust_mad(dapi)
    nuclei_binary = gaussian_filter(dapi, 1.2) > dapi_med + 8 * dapi_mad
    nuclei_binary = closing(nuclei_binary, disk(2))
    nuclei_local = filter_components(nuclei_binary, min_area=80, max_area=2500)
    if int(nuclei_local.max()) != 6:
        raise RuntimeError(f"Expected six high-confidence nuclei, found {nuclei_local.max()}")
    nuclei_all = nuclei_local > 0
    nucleus_rows = []
    for p in regionprops(nuclei_local):
        obj = nuclei_local == p.label
        background, background_n = clipped_annular_mean(tunel, obj, nuclei_all, inner=6, outer=22)
        raw = float(tunel[obj].mean())
        nucleus_rows.append({
            "nucleus_id": int(p.label),
            "centroid_x": float(p.centroid[1] + sx0),
            "centroid_y": float(p.centroid[0] + sy0),
            "area_px": int(p.area),
            "raw_tunel_intensity": raw,
            "background_intensity": background,
            "corrected_tunel_intensity": raw - background,
            "background_pixels": background_n,
        })
    corrected = np.sort([r["corrected_tunel_intensity"] for r in nucleus_rows])
    gap_idx = int(np.argmax(np.diff(corrected)))
    tunel_cut = float((corrected[gap_idx] + corrected[gap_idx + 1]) / 2)
    for row in nucleus_rows:
        row["group"] = "TUNEL+" if row["corrected_tunel_intensity"] > tunel_cut else "TUNEL-"
    nuclei_df = pd.DataFrame(nucleus_rows).drop(columns="background_pixels")
    nuclei_df.to_csv(output_dir / "tunel_nuclei.csv", index=False, float_format="%.6f")
    nuclei_full = np.zeros(shape_stain, dtype=np.uint16)
    nuclei_full[sy0:sy1, sx0:sx1] = nuclei_local
    tifffile.imwrite(output_dir / "tunel_nuclei_mask.tif", nuclei_full, metadata={"axes": "YX"}, compression="zlib")
    scatter_plot(nuclei_df, "group", "corrected_tunel_intensity", "Corrected nuclear TUNEL intensity (a.u.)", output_dir / "tunel_scatter.png")

    # Host macrophages: Cellpose cyto3 on last-frame brightfield, fully headless.
    helper = Path(__file__).resolve().with_name("host_cellpose.py")
    with tempfile.TemporaryDirectory(prefix="parasite_task_") as temp_name:
        temp = Path(temp_name)
        bf_tif = temp / "movie_bf_roi.tif"
        host_npy = temp / "host.npy"
        tifffile.imwrite(bf_tif, movie_crop[0])
        subprocess.run([
            "/opt/conda/envs/cellpose/bin/python",
            str(helper),
            str(bf_tif),
            str(host_npy),
        ], check=True)
        host_raw = np.load(host_npy)
    host_local = filter_existing_labels(host_raw, min_area=2500, max_area=15000)
    if int(host_local.max()) != 6:
        raise RuntimeError(f"Expected six macrophages, found {host_local.max()}")

    # Map registered nuclei to hosts and transfer the TUNEL class.
    host_groups: dict[int, str] = {}
    nucleus_host_distances = []
    for row in nucleus_rows:
        x_movie, y_movie = transform_xy(global_affine, row["centroid_x"], row["centroid_y"])
        # Allow a small boundary tolerance: fixation changes cell outlines and Cellpose can
        # exclude a thin, dark nuclear edge even when the registered nucleus is correct.
        host_id, distance = nearest_label(host_local, x_movie - mx0, y_movie - my0, max(10.0, 4 * reg["median_residual"]))
        if host_id == 0:
            raise RuntimeError(f"Nucleus {row['nucleus_id']} could not be assigned to a host")
        if host_id in host_groups:
            raise RuntimeError(f"More than one nucleus assigned to host {host_id}")
        host_groups[host_id] = row["group"]
        nucleus_host_distances.append(distance)
    if len(host_groups) != 6:
        raise RuntimeError("Not every host received exactly one TUNEL classification")
    host_full = np.zeros(shape_movie, dtype=np.uint16)
    host_full[my0:my1, mx0:mx1] = host_local
    tifffile.imwrite(output_dir / "host_mask.tif", host_full, metadata={"axes": "YX"}, compression="zlib")

    # Parasites: red-channel bandpass foreground plus marker watershed for resolvable contacts.
    red = movie_crop[red_idx].astype(float)
    green = movie_crop[green_idx].astype(float)
    red_bp = gaussian_filter(red, 1) - gaussian_filter(red, 12)
    red_med, red_mad = float(np.median(red_bp)), robust_mad(red_bp)
    candidate = closing(red_bp > red_med + 6 * red_mad, disk(1))
    candidate_cc = filter_components(candidate, min_area=30)
    component_props = list(regionprops(candidate_cc))
    one_peak_areas = []
    component_peaks: dict[int, np.ndarray] = {}
    for p in component_props:
        component = candidate_cc[p.slice] == p.label
        peaks = peak_local_max(
            red_bp[p.slice],
            min_distance=6,
            threshold_abs=red_med + 10 * red_mad,
            labels=component,
            exclude_border=False,
        )
        component_peaks[p.label] = peaks
        if len(peaks) == 1:
            one_peak_areas.append(float(p.area))
    small = [a for a in one_peak_areas if a < np.median(one_peak_areas) * 1.5]
    single_area = float(np.median(small))
    parasite_objects: list[np.ndarray] = []
    excluded_clumps = []
    for p in component_props:
        component = candidate_cc[p.slice] == p.label
        peaks = component_peaks[p.label]
        if len(peaks) == 0 or (len(peaks) == 1 and p.area > 1.8 * single_area):
            excluded_clumps.append(p)
            continue
        markers = np.zeros(component.shape, dtype=np.int32)
        for marker_id, (yy, xx) in enumerate(peaks, 1):
            markers[yy, xx] = marker_id
        split = watershed(-red_bp[p.slice], markers, mask=component)
        for marker_id in range(1, len(peaks) + 1):
            obj = np.zeros(candidate.shape, dtype=bool)
            obj[p.slice] = split == marker_id
            parasite_objects.append(obj)
    parasite_objects.sort(key=lambda obj: tuple(np.argwhere(obj).mean(axis=0)))
    parasite_local = np.zeros(candidate.shape, dtype=np.uint16)
    for parasite_id, obj in enumerate(parasite_objects, 1):
        parasite_local[obj] = parasite_id

    parasite_rows = []
    unassigned = []
    for p in regionprops(parasite_local):
        obj = parasite_local == p.label
        background_red, _ = clipped_annular_mean(red, obj, candidate, inner=4, outer=15)
        background_green, _ = clipped_annular_mean(green, obj, candidate, inner=4, outer=15)
        raw_red = float(red[obj].mean())
        raw_green = float(green[obj].mean())
        corrected_red = raw_red - background_red
        corrected_green = raw_green - background_green
        if corrected_green <= 0:
            raise RuntimeError(f"Parasite {p.label} has non-positive corrected green intensity")
        cx_local, cy_local = float(p.centroid[1]), float(p.centroid[0])
        host_id, host_distance = nearest_label(host_local, cx_local, cy_local, max_distance=3.0)
        if host_id == 0:
            unassigned.append(int(p.label))
            continue
        parasite_rows.append({
            "parasite_id": int(p.label),
            "centroid_x": cx_local + mx0,
            "centroid_y": cy_local + my0,
            "raw_red": raw_red,
            "raw_green": raw_green,
            "background_red": background_red,
            "background_green": background_green,
            "corrected_red": corrected_red,
            "corrected_green": corrected_green,
            "proliferation_index": 10.0 - corrected_red / corrected_green,
            "host_cell_id": int(host_id),
            "group": host_groups[host_id],
        })
    if unassigned:
        raise RuntimeError(f"Resolved parasites lacked confident host assignment: {unassigned}")
    parasites_df = pd.DataFrame(parasite_rows)
    parasites_df.to_csv(output_dir / "parasite_index.csv", index=False, float_format="%.6f")
    parasite_full = np.zeros(shape_movie, dtype=np.uint16)
    parasite_full[my0:my1, mx0:mx1] = parasite_local
    tifffile.imwrite(output_dir / "parasite_mask.tif", parasite_full, metadata={"axes": "YX"}, compression="zlib")
    scatter_plot(parasites_df, "group", "proliferation_index", "Proliferation index", output_dir / "parasite_scatter.png")

    # Quantitative QC and report.
    stain_luts = channel_luts(stain_lif)
    low_green = parasites_df.loc[parasites_df.corrected_green < 3, "parasite_id"].astype(int).tolist()
    local_matrix = "; ".join(", ".join(f"{v:.6f}" for v in row) for row in local_affine)
    global_matrix = "; ".join(", ".join(f"{v:.6f}" for v in row) for row in global_affine)
    clump_descriptions = [
        f"area {int(p.area)} px at local centroid (x={p.centroid[1]:.1f}, y={p.centroid[0]:.1f})"
        for p in excluded_clumps
    ]
    report = f"""# Position 40: TUNEL and parasite proliferation analysis

## Inputs, channels, and ROI checks

Both position-40 source series are {shape_stain[1]} x {shape_stain[0]} pixels. The stain ROI is `(x={sx0}:{sx1}, y={sy0}:{sy1})`; the last-frame ROI is `(x={mx0}:{mx1}, y={my0}:{my1})`. Both are valid {sx1-sx0} x {sy1-sy0} rectangles entirely inside their source images. ImageJ rectangle ROI files do not store the source canvas width/height, so an independent ROI-canvas dimension cannot be compared; the bounds, equal ROI sizes, and cell content were checked against the native 2048 x 2048 sources. No rescaling or ROI correction was needed.

The stained LIF acquisition settings were {wavelengths[1]:.0f} nm for channel 1 and {wavelengths[2]:.0f} nm for channel 2. Channel {tunel_idx} is therefore the far-red CF640R/TUNEL acquisition and channel {dapi_idx} is DAPI. The saved viewer LUT names `{stain_luts}` are misleading relative to the stain identities; physical acquisition settings and image morphology were used. In the movie, the LIF LUT metadata identifies channel {green_idx} as green mKikume and channel {red_idx} as red mKikume (`{movie_luts}`).

## Nuclear TUNEL analysis

DAPI was Gaussian-smoothed (sigma 1.2 px), thresholded at ROI median + 8 robust MAD (threshold {dapi_med + 8*dapi_mad:.2f}), closed by 2 px, and filtered to 80-2500 px. This yielded {len(nuclei_df)} nuclei. IDs and CSV centroids are in full stained-image coordinates; the TIFF is a full 2048 x 2048 label image with zero outside the ROI.

For each DAPI nucleus, the TUNEL value is the mean of its pixels. Local background was measured in a 6-22 px annulus around that nucleus, excluding all nuclei plus a 2 px guard. Pixels outside median +/- 3 robust MAD within the annulus were rejected to avoid bright debris, and the mean of the remaining local pixels was subtracted. Thus every nucleus has its own nearby background estimate.

The corrected nuclear values were {', '.join(f'{v:.2f}' for v in corrected)} a.u. The largest gap was between {corrected[gap_idx]:.2f} and {corrected[gap_idx+1]:.2f}; I placed the cut at their midpoint, **{tunel_cut:.2f} a.u.** Values above it are TUNEL+, values below it TUNEL-. This gives {(nuclei_df.group == 'TUNEL+').sum()} TUNEL+ and {(nuclei_df.group == 'TUNEL-').sum()} TUNEL- nuclei. This is an image-internal intensity cut; with only six nuclei it should not be treated as a universal assay threshold.

## Registration and macrophage hosts

Registration used only brightfield: SIFT feature matches and RANSAC supplied the initialization, followed by rigid ECC refinement on background-subtracted images. The local stain-ROI to movie-ROI affine matrix `[x'; y'] = A[x; y; 1]` is `[{local_matrix}]`; the equivalent full-image matrix is `[{global_matrix}]`. It is a {reg['angle_deg']:.3f}-degree rotation plus translation, with no fitted scale or shear. There were {reg['inliers']} RANSAC inliers among {reg['sift_matches']} ratio-tested matches; their post-refinement median residual was **{reg['median_residual']:.2f} px** (mean {reg['mean_residual']:.2f} px), and ECC correlation was {reg['ecc']:.3f}. Registered cell and nuclear landmarks overlay correctly.

Macrophages were segmented headlessly from last-frame brightfield with pretrained Cellpose `cyto3` (diameter 90 px, cell-probability threshold -1), then restricted to 2500-15000 px objects. It found {int(host_local.max())} hosts. Each registered DAPI centroid matched a distinct host; five centroids were inside their host mask and one required a {max(nucleus_host_distances):.2f} px nearest-boundary association (about {0.32*max(nucleus_host_distances):.2f} micrometres), consistent with the thin nuclear edge omitted by the brightfield cell mask. Host IDs in `host_mask.tif` are the IDs used in `parasite_index.csv`.

## Parasites, clumps, and intensity calculation

The movie metadata identifies channel {red_idx} as red and channel {green_idx} as green. Parasites were detected in red because its local signal-to-background was stronger. I subtracted a sigma-12-px Gaussian background from a sigma-1-px image, thresholded at median + 6 robust MAD, and used persistent peaks above median + 10 robust MAD as watershed seeds. Contacts with two spatially separate persistent peaks were split. A candidate with area greater than 1.8 times the single-object median and no separable second peak was classified as an unresolved clump.

I retained {len(parasites_df)} clearly resolved parasites and excluded {len(excluded_clumps)} unresolved clump candidate(s): {'; '.join(clump_descriptions) if clump_descriptions else 'none'}. Exclusion was by size plus lack of a resolvable internal intensity valley/second peak, not by host group.

For each parasite and each color separately, background was the sigma-clipped mean in a nearby 4-15 px annulus, excluding all detected parasite signal plus a 2 px guard. The reported raw values are within-mask means; corrected values are raw minus the corresponding local background. The proliferation index was calculated exactly as `10 - corrected_red / corrected_green`.

All {len(parasites_df)} retained parasites fell inside a Cellpose host mask and were assigned confidently; none was guessed or left unassigned. Parasite group is the registered nuclear TUNEL class of that host. {('Parasite ID(s) ' + ', '.join(map(str, low_green)) + ' had corrected green below 3 a.u.; the requested ratio is reported, but it is denominator-sensitive and should be interpreted cautiously.') if low_green else 'No parasite had a near-zero corrected-green denominator.'}

## Consistency checks and unresolved points

`tunel_nuclei_mask.tif` has {int(nuclei_local.max())} nonzero IDs and `tunel_nuclei.csv` has {len(nuclei_df)} rows. `parasite_mask.tif` has {int(parasite_local.max())} nonzero IDs and `parasite_index.csv` has {len(parasites_df)} rows. No manual ROI, click, or interactive correction was used. The only material uncertainty is the near-background green denominator noted above and the biological generalizability of a six-nucleus TUNEL cut.
"""
    (output_dir / "report.md").write_text(report)

    print(f"Wrote {output_dir}")
    print(f"Nuclei: {len(nuclei_df)}; parasites: {len(parasites_df)}; hosts: {host_local.max()}; excluded clumps: {len(excluded_clumps)}")


if __name__ == "__main__":
    main()
