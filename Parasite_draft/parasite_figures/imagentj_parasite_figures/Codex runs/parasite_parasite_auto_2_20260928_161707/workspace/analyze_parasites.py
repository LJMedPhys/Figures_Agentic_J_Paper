#!/usr/bin/env python3
"""Fully automatic position-40 single-parasite proliferation analysis."""

from __future__ import annotations

import math
import struct
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from readlif.reader import LifFile
from scipy.ndimage import distance_transform_edt, gaussian_filter, label as ndi_label
from scipy.optimize import linear_sum_assignment
from sklearn.mixture import GaussianMixture
from skimage.feature import peak_local_max
from skimage.measure import regionprops
from skimage.morphology import dilation, disk
from skimage.segmentation import find_boundaries, watershed


ROOT = Path("data/parasite_task")
OUT = Path("output")
DIAG = OUT / "_diagnostics"
MOVIE_PATH = ROOT / "time-lapse/20211028_Pos36_40.lif"
FIXED_PATH = ROOT / "20210728_TUNEL_staining.lif"
ROI_PATH = ROOT / "one_frame/LastFrame_Overview.roi"
CELLPOSE_PYTHON = Path("/opt/conda/envs/cellpose/bin/python")
CELLPOSE_SCRIPT = Path(__file__).with_name("segment_hosts_cellpose.py")


@dataclass
class Nucleus:
    label: int
    centroid_fixed_xy: tuple[float, float]
    centroid_roi_xy: tuple[float, float]
    tunel_corrected: float
    tunel_group: str = ""


def read_imagej_rect(path: Path) -> tuple[int, int, int, int]:
    data = path.read_bytes()
    if data[:4] != b"Iout":
        raise ValueError(f"Not an ImageJ ROI: {path}")
    roi_type = data[6]
    if roi_type != 1:
        raise ValueError("The supplied analysis ROI must be rectangular")
    top, left, bottom, right = struct.unpack(">4H", data[8:16])
    return left, top, right, bottom


def position_channel_metadata(lif: LifFile, position: str) -> list[dict]:
    root = ET.fromstring(lif.xml_header)
    target = None
    for element in root.iter():
        if element.attrib.get("Name") == position:
            target = element
            break
    if target is None:
        raise RuntimeError(f"Could not find metadata for {position}")
    channels = []
    for element in target.iter():
        if element.tag.split("}")[-1] == "WideFieldChannelInfo":
            a = element.attrib
            if a.get("Active") == "1":
                channels.append({
                    "lut": a.get("LUT", ""),
                    "emission_nm": float(a.get("EmissionWavelength", "nan")),
                    "method": a.get("ContrastingMethodName", ""),
                    "name": a.get("UserDefName", ""),
                })
    if len(channels) < 3:
        raise RuntimeError("Incomplete channel metadata")
    return channels[:3]


def robust_u8(image: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(image, [1, 99])
    return np.uint8(255 * np.clip((image.astype(float) - lo) / max(hi - lo, 1), 0, 1))


def robust01(image: np.ndarray, lo=1, hi=99) -> np.ndarray:
    a, b = np.percentile(image, [lo, hi])
    return np.clip((image.astype(float) - a) / max(b - a, 1), 0, 1)


def register_brightfield(fixed_bf: np.ndarray, movie_bf: np.ndarray):
    """Return a fixed->movie similarity transform and landmark diagnostics."""
    cv2.setRNGSeed(20210728)
    target, moving = robust_u8(movie_bf), robust_u8(fixed_bf)
    sift = cv2.SIFT_create(nfeatures=5000, contrastThreshold=0.02, edgeThreshold=15)
    kp_target, des_target = sift.detectAndCompute(target, None)
    kp_moving, des_moving = sift.detectAndCompute(moving, None)
    pairs = cv2.BFMatcher().knnMatch(des_moving, des_target, k=2)
    matches = [m for m, n in pairs if m.distance < 0.72 * n.distance]
    src = np.float32([kp_moving[m.queryIdx].pt for m in matches])
    dst = np.float32([kp_target[m.trainIdx].pt for m in matches])
    matrix, inliers = cv2.estimateAffinePartial2D(
        src, dst, method=cv2.RANSAC, ransacReprojThreshold=4,
        maxIters=10000, confidence=0.999,
    )
    if matrix is None or int(inliers.sum()) < 20:
        raise RuntimeError("Brightfield registration failed quality control")
    predicted = cv2.transform(src[None], matrix)[0]
    residuals = np.linalg.norm(predicted - dst, axis=1)
    good = inliers.ravel().astype(bool)
    return matrix, len(matches), int(good.sum()), residuals[good]


def dog_snr(image: np.ndarray) -> tuple[np.ndarray, float]:
    dog = gaussian_filter(image.astype(float), 1) - gaussian_filter(image.astype(float), 12)
    median = np.median(dog)
    sigma = 1.4826 * np.median(np.abs(dog - median))
    return dog / max(sigma, 1e-6), sigma


def segment_parasites(green: np.ndarray, red: np.ndarray):
    """Detect clear singles, split well-separated objects, reject unresolved clumps."""
    zg, noise_g = dog_snr(green)
    zr, noise_r = dog_snr(red)
    score = np.maximum(zg, zr)
    primary = score > 5.0
    components, _ = ndi_label(primary)
    accepted = []
    clumps = []
    boundary_rejections = []
    target_objects = []

    for region in regionprops(components, intensity_image=score):
        y0, x0, y1, x1 = region.bbox
        # A two-pixel guard is needed because the Gaussian PSF can terminate just
        # inside a cropped boundary even when the biological object is truncated.
        touches_border = y0 <= 2 or x0 <= 2 or y1 >= score.shape[0] - 2 or x1 >= score.shape[1] - 2
        if touches_border and region.intensity_max >= 8:
            boundary_rejections.append({"bbox": region.bbox, "area": int(region.area)})
            continue
        if region.area < 15:
            continue
        component = components == region.label
        target_objects.append(component)
        if region.area <= 200:
            accepted.append(component)
            continue

        peaks = peak_local_max(
            score, min_distance=3, threshold_abs=8,
            labels=component, exclude_border=False,
        )
        if len(peaks) < 2:
            clumps.append({"bbox": region.bbox, "area": int(region.area), "peaks": len(peaks)})
            continue
        distances = np.linalg.norm(peaks[:, None, :] - peaks[None, :, :], axis=2)
        distances[distances == 0] = np.inf
        minimum_separation = float(distances.min())
        if minimum_separation < 12:
            clumps.append({
                "bbox": region.bbox, "area": int(region.area),
                "peaks": len(peaks), "minimum_peak_separation": minimum_separation,
            })
            continue

        markers = np.zeros_like(components, dtype=np.int32)
        for marker_id, (py, px) in enumerate(peaks, 1):
            markers[py, px] = marker_id
        split = watershed(-score, markers, mask=component)
        for split_id in range(1, int(split.max()) + 1):
            part = split == split_id
            if part.sum() >= 15:
                accepted.append(part)

    accepted.sort(key=lambda m: tuple(np.argwhere(m).mean(axis=0)))
    labels = np.zeros(score.shape, dtype=np.uint16)
    for parasite_id, mask in enumerate(accepted, 1):
        labels[mask] = parasite_id
    return labels, score, zg, zr, noise_g, noise_r, clumps, boundary_rejections, target_objects


def segment_and_measure_nuclei(
    dapi: np.ndarray, tunel: np.ndarray, transform: np.ndarray,
    roi: tuple[int, int, int, int],
):
    """Segment all nuclei and use their local TUNEL signal to learn two groups."""
    top_hat = gaussian_filter(dapi.astype(float), 1) - gaussian_filter(dapi.astype(float), 20)
    binary = top_hat > 100
    initial, _ = ndi_label(binary)
    valid = np.zeros_like(initial, dtype=np.int32)
    valid_regions = []
    for region in regionprops(initial):
        if 30 <= region.area <= 1500 and region.solidity > 0.5:
            valid_regions.append(region)
            valid[initial == region.label] = len(valid_regions)

    full_records = []
    for new_id, region in enumerate(valid_regions, 1):
        mask = valid == new_id
        cy, cx = region.centroid
        pad = 24
        y0, y1 = max(0, int(cy) - pad), min(dapi.shape[0], int(cy) + pad + 1)
        x0, x1 = max(0, int(cx) - pad), min(dapi.shape[1], int(cx) + pad + 1)
        yy, xx = np.ogrid[y0:y1, x0:x1]
        local_mask = mask[y0:y1, x0:x1]
        annulus = ((yy - cy) ** 2 + (xx - cx) ** 2 >= 10 ** 2)
        annulus &= ((yy - cy) ** 2 + (xx - cx) ** 2 <= 22 ** 2)
        annulus &= valid[y0:y1, x0:x1] == 0
        corrected = float(tunel[mask].mean() - np.median(tunel[y0:y1, x0:x1][annulus]))
        full_records.append((new_id, float(cx), float(cy), corrected))

    values = np.array([r[3] for r in full_records])
    mixture = GaussianMixture(2, random_state=20210728, n_init=20).fit(values[:, None])
    means = mixture.means_.ravel()
    low_component, high_component = np.argsort(means)
    grid = np.linspace(means[low_component], means[high_component], 2000)
    posterior = mixture.predict_proba(grid[:, None])
    threshold = float(grid[np.argmin(np.abs(posterior[:, low_component] - posterior[:, high_component]))])

    x0, y0, x1, y1 = roi
    roi_nuclei = []
    for nucleus_id, cx, cy, corrected in full_records:
        mapped = cv2.transform(np.float32([[[cx, cy]]]), transform)[0, 0]
        rx, ry = float(mapped[0] - x0), float(mapped[1] - y0)
        if 0 <= rx < x1 - x0 and 0 <= ry < y1 - y0:
            group = "TUNEL+" if mixture.predict([[corrected]])[0] == high_component else "TUNEL-"
            roi_nuclei.append(Nucleus(nucleus_id, (cx, cy), (rx, ry), corrected, group))
    return roi_nuclei, threshold, means[np.argsort(means)], len(full_records)


def closest_candidate_label(mask: np.ndarray, xy: tuple[float, float]) -> tuple[int, float]:
    x, y = xy
    iy = int(np.clip(round(y), 0, mask.shape[0] - 1))
    ix = int(np.clip(round(x), 0, mask.shape[1] - 1))
    direct = int(mask[iy, ix])
    if direct:
        return direct, 0.0
    best_label, best_distance = 0, np.inf
    for region in regionprops(mask):
        y0, x0, y1, x1 = region.bbox
        if x < x0 - 40 or x > x1 + 40 or y < y0 - 40 or y > y1 + 40:
            continue
        local = mask[max(0, y0 - 40):min(mask.shape[0], y1 + 40),
                     max(0, x0 - 40):min(mask.shape[1], x1 + 40)] == region.label
        oy, ox = max(0, y0 - 40), max(0, x0 - 40)
        distance = distance_transform_edt(~local)[iy - oy, ix - ox]
        if distance < best_distance:
            best_label, best_distance = region.label, float(distance)
    return int(best_label), float(best_distance)


def select_and_match_hosts(
    candidates: np.ndarray, parasite_mask: np.ndarray, target_objects: list[np.ndarray],
    nuclei: list[Nucleus],
):
    """Filter Cellpose candidates using infection signal and match them to nuclei."""
    selected = set()
    infection_objects = [parasite_mask == i for i in range(1, int(parasite_mask.max()) + 1)]
    infection_objects += target_objects
    for obj in infection_objects:
        if not obj.any():
            continue
        cy, cx = np.argwhere(obj).mean(axis=0)
        candidate, distance = closest_candidate_label(candidates, (cx, cy))
        if candidate and distance <= 30:
            selected.add(candidate)

    candidate_regions = {r.label: r for r in regionprops(candidates)}
    nucleus_xy = np.array([n.centroid_roi_xy[::-1] for n in nuclei])  # y, x
    selected_list = sorted(selected)
    if not selected_list:
        raise RuntimeError("No Cellpose host candidates contain infection signal")
    selected_xy = np.array([candidate_regions[i].centroid for i in selected_list])
    cost = np.linalg.norm(selected_xy[:, None, :] - nucleus_xy[None, :, :], axis=2)
    row, col = linear_sum_assignment(cost)
    matches = {selected_list[r]: int(c) for r, c in zip(row, col)}

    remaining_nuclei = [i for i in range(len(nuclei)) if i not in matches.values()]
    remaining_candidates = [i for i in candidate_regions if i not in matches]
    for nucleus_index in remaining_nuclei:
        if not remaining_candidates:
            break
        distances = [
            np.linalg.norm(np.array(candidate_regions[i].centroid) - nucleus_xy[nucleus_index])
            for i in remaining_candidates
        ]
        best = int(np.argmin(distances))
        if distances[best] <= 90:
            label = remaining_candidates.pop(best)
            matches[label] = nucleus_index

    if len(matches) != len(nuclei):
        raise RuntimeError(f"Only matched {len(matches)} of {len(nuclei)} macrophages")

    ordered = sorted(matches.items(), key=lambda item: nuclei[item[1]].centroid_roi_xy[1])
    host_mask = np.zeros_like(candidates, dtype=np.uint16)
    host_records = []
    for host_id, (candidate_label, nucleus_index) in enumerate(ordered, 1):
        host_mask[candidates == candidate_label] = host_id
        nucleus = nuclei[nucleus_index]
        host_records.append({
            "host_id": host_id,
            "candidate_label": candidate_label,
            "nucleus": nucleus,
            "centroid_yx": candidate_regions[candidate_label].centroid,
        })
    return host_mask, host_records


def parasite_measurements(
    parasite_mask: np.ndarray, host_mask: np.ndarray, red: np.ndarray, green: np.ndarray,
    score: np.ndarray, host_records: list[dict], roi_origin: tuple[int, int],
):
    rows, uncertain = [], []
    group_by_host = {r["host_id"]: r["nucleus"].tunel_group for r in host_records}
    all_fluorescence = dilation(score >= 3, disk(3))
    for parasite_id in range(1, int(parasite_mask.max()) + 1):
        mask = parasite_mask == parasite_id
        cy, cx = np.argwhere(mask).mean(axis=0)
        raw_red, raw_green = float(red[mask].mean()), float(green[mask].mean())
        annulus = dilation(mask, disk(15)) & ~dilation(mask, disk(6))
        annulus &= ~all_fluorescence
        annulus &= parasite_mask == 0
        if annulus.sum() < 100:
            raise RuntimeError(f"Insufficient local background for parasite {parasite_id}")
        background_red = float(red[annulus].mean())
        background_green = float(green[annulus].mean())
        corrected_red = raw_red - background_red
        corrected_green = raw_green - background_green
        if corrected_green <= 0:
            raise RuntimeError(f"Nonpositive corrected green signal for parasite {parasite_id}")

        overlaps, counts = np.unique(host_mask[mask], return_counts=True)
        nonzero = overlaps > 0
        if nonzero.any():
            host_id = int(overlaps[nonzero][np.argmax(counts[nonzero])])
            overlap_fraction = counts[nonzero].max() / mask.sum()
            if overlap_fraction < 0.5:
                uncertain.append(parasite_id)
        else:
            host_id = 0
            uncertain.append(parasite_id)
        rows.append({
            "parasite_id": parasite_id,
            "centroid_x": cx + roi_origin[0],
            "centroid_y": cy + roi_origin[1],
            "raw_red": raw_red,
            "raw_green": raw_green,
            "background_red": background_red,
            "background_green": background_green,
            "corrected_red": corrected_red,
            "corrected_green": corrected_green,
            "proliferation_index": 10 - corrected_red / corrected_green,
            "host_cell_id": host_id if host_id else pd.NA,
            "group": group_by_host.get(host_id, pd.NA),
        })
    return pd.DataFrame(rows), uncertain


def save_diagnostics(
    last_roi, fixed_roi, parasite_mask, host_mask, score, clumps, dataframe,
):
    # Registration edge overlay.
    def edges(x):
        x = gaussian_filter(x.astype(float), 1)
        gx = cv2.Sobel(x, cv2.CV_64F, 1, 0)
        gy = cv2.Sobel(x, cv2.CV_64F, 0, 1)
        e = np.hypot(gx, gy)
        return e / max(np.percentile(e, 99), 1)

    live_edge, fixed_edge = edges(last_roi[0]), edges(fixed_roi[0])
    overlay = np.stack([live_edge, fixed_edge, np.zeros_like(live_edge)], axis=-1)
    plt.imsave(DIAG / "registration_edges.png", np.clip(overlay, 0, 1))

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    axes[0].imshow(robust01(last_roi[0]), cmap="gray")
    axes[0].imshow(find_boundaries(host_mask), cmap="Reds", alpha=0.8)
    axes[0].set_title("last BF: retained Cellpose host boundaries")
    axes[1].imshow(score, cmap="magma", vmin=0, vmax=min(25, np.percentile(score, 99.9)))
    axes[1].imshow(find_boundaries(parasite_mask), cmap="winter", alpha=0.9)
    axes[1].set_title("combined fluorescence SNR and parasite masks")
    axes[2].imshow(robust01(fixed_roi[2], 2, 99.8), cmap="Blues")
    axes[2].imshow(robust01(fixed_roi[1], 2, 99.8), cmap="Greens", alpha=0.55)
    axes[2].set_title("registered DAPI (blue) + TUNEL (green)")
    for row in dataframe.itertuples():
        x, y = row.centroid_x - 537, row.centroid_y - 849
        axes[1].text(x, y, str(row.parasite_id), color="cyan", fontsize=9)
    for clump in clumps:
        y0, x0, y1, x1 = clump["bbox"]
        axes[1].add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                        fill=False, edgecolor="yellow", linewidth=2))
    for record in range(1, int(host_mask.max()) + 1):
        cy, cx = np.argwhere(host_mask == record).mean(axis=0)
        axes[0].text(cx, cy, f"H{record}", color="yellow", fontsize=10)
    for ax in axes:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(DIAG / "segmentation_qc.png", dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    DIAG.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1 = read_imagej_rect(ROI_PATH)
    roi = (x0, y0, x1, y1)

    movie_lif, fixed_lif = LifFile(str(MOVIE_PATH)), LifFile(str(FIXED_PATH))
    movie_meta = position_channel_metadata(movie_lif, "Position040")
    fixed_meta = position_channel_metadata(fixed_lif, "Position040")
    movie = movie_lif.get_image(4)
    fixed = fixed_lif.get_image(39)
    last = np.stack([np.asarray(movie.get_frame(t=movie.nt - 1, c=c)) for c in range(3)])
    stained = np.stack([np.asarray(fixed.get_frame(c=c)) for c in range(3)])

    movie_bf = next(i for i, m in enumerate(movie_meta) if m["method"] == "TL-BF")
    movie_fluorescence = [i for i in range(3) if i != movie_bf]
    green_channel = min(movie_fluorescence, key=lambda i: abs(movie_meta[i]["emission_nm"] - 525))
    red_channel = min(movie_fluorescence, key=lambda i: abs(movie_meta[i]["emission_nm"] - 600))
    fixed_bf = next(i for i, m in enumerate(fixed_meta) if m["method"] == "TL-BF")
    fixed_fluorescence = [i for i in range(3) if i != fixed_bf]
    tunel_channel = max(fixed_fluorescence, key=lambda i: fixed_meta[i]["emission_nm"])
    dapi_channel = min(fixed_fluorescence, key=lambda i: fixed_meta[i]["emission_nm"])

    transform, total_matches, inliers, residuals = register_brightfield(
        stained[fixed_bf], last[movie_bf],
    )
    registered = np.stack([
        cv2.warpAffine(channel, transform, (last.shape[2], last.shape[1]), flags=cv2.INTER_LINEAR)
        for channel in stained
    ])
    last_roi = last[:, y0:y1, x0:x1]
    fixed_roi = registered[:, y0:y1, x0:x1]
    tifffile.imwrite(DIAG / "last_roi.tif", last_roi, metadata={"axes": "CYX"})
    tifffile.imwrite(DIAG / "fixed_registered_roi.tif", fixed_roi, metadata={"axes": "CYX"})

    parasite_mask, score, zg, zr, noise_g, noise_r, clumps, boundary_rejections, target_objects = segment_parasites(
        last_roi[green_channel], last_roi[red_channel],
    )

    cellpose_mask_path = DIAG / "cellpose_host_candidates.tif"
    subprocess.run([
        str(CELLPOSE_PYTHON), str(CELLPOSE_SCRIPT),
        str(DIAG / "last_roi.tif"), str(cellpose_mask_path),
    ], check=True)
    host_candidates = tifffile.imread(cellpose_mask_path)

    nuclei, tunel_threshold, tunel_means, total_nuclei = segment_and_measure_nuclei(
        stained[dapi_channel], stained[tunel_channel], transform, roi,
    )
    host_mask_roi, host_records = select_and_match_hosts(
        host_candidates, parasite_mask, target_objects, nuclei,
    )
    dataframe, uncertain = parasite_measurements(
        parasite_mask, host_mask_roi, last_roi[red_channel], last_roi[green_channel],
        score, host_records, (x0, y0),
    )

    # Full-frame label masks keep CSV coordinates in the original position-40 coordinate system.
    parasite_full = np.zeros(last.shape[1:], dtype=np.uint16)
    host_full = np.zeros(last.shape[1:], dtype=np.uint16)
    parasite_full[y0:y1, x0:x1] = parasite_mask
    host_full[y0:y1, x0:x1] = host_mask_roi
    tifffile.imwrite(OUT / "parasite_mask.tif", parasite_full)
    tifffile.imwrite(OUT / "host_mask.tif", host_full)
    numeric = [c for c in dataframe.columns if c not in {"parasite_id", "host_cell_id", "group"}]
    dataframe[numeric] = dataframe[numeric].round(6)
    dataframe.to_csv(OUT / "parasite_index.csv", index=False)

    # Requested grouped dot plot.
    groups = ["TUNEL-", "TUNEL+"]
    colors = {"TUNEL-": "#3274A1", "TUNEL+": "#D64B40"}
    rng = np.random.default_rng(20210728)
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    for gx, group in enumerate(groups):
        values = dataframe.loc[dataframe.group == group, "proliferation_index"].to_numpy()
        ax.scatter(gx + rng.uniform(-0.09, 0.09, len(values)), values, s=38,
                   color=colors[group], edgecolor="white", linewidth=0.6, zorder=3)
        if len(values):
            median = float(np.median(values))
            ax.plot([gx - 0.22, gx + 0.22], [median, median], color="black", lw=2.2, zorder=4)
    ax.set_xticks([0, 1], groups)
    ax.set_ylabel("Proliferation index")
    ax.set_xlim(-0.5, 1.5)
    ax.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(OUT / "parasite_scatter.png", dpi=220)
    plt.close(fig)

    save_diagnostics(last_roi, fixed_roi, parasite_mask, host_mask_roi, score, clumps, dataframe)

    scale = math.sqrt(transform[0, 0] ** 2 + transform[1, 0] ** 2)
    angle = math.degrees(math.atan2(transform[1, 0], transform[0, 0]))
    host_lines = []
    for record in host_records:
        n = record["nucleus"]
        count = int((dataframe.host_cell_id == record["host_id"]).sum())
        host_lines.append(
            f"| {record['host_id']} | {n.tunel_corrected:.2f} | {n.tunel_group} | {count} |"
        )
    median_lines = []
    for group in groups:
        vals = dataframe.loc[dataframe.group == group, "proliferation_index"]
        median_lines.append(f"- {group}: n={len(vals)}, median={vals.median():.3f}")
    uncertain_text = (
        "None; every retained parasite overlapped one retained host mask by at least 50% of its pixels."
        if not uncertain else f"Parasite IDs {uncertain} had <50% host-mask overlap and remain unassigned/flagged."
    )
    clump_details = "; ".join(
        f"bbox (x={c['bbox'][1]}:{c['bbox'][3]-1}, y={c['bbox'][0]}:{c['bbox'][2]-1}), "
        f"area={c['area']} px, peaks={c['peaks']}, minimum peak spacing={c.get('minimum_peak_separation', float('nan')):.1f} px"
        for c in clumps
    ) or "none"
    boundary_text = (
        f"In addition, {len(boundary_rejections)} strong fluorescence object at the left ROI boundary "
        "was rejected as boundary-truncated (not counted as a clump), because its full shape could not be evaluated."
        if len(boundary_rejections) == 1 else
        (f"In addition, {len(boundary_rejections)} strong fluorescence objects at ROI boundaries were rejected "
         "as boundary-truncated (not counted as clumps), because their full shapes could not be evaluated."
         if boundary_rejections else
         "No strong candidate touched an ROI boundary, so no boundary-truncated object required exclusion.")
    )
    report = f"""# Position 40, last-frame parasite proliferation analysis

## Scope and coordinate convention

The supplied `LastFrame_Overview.ROI` was read directly as x={x0}:{x1 - 1}, y={y0}:{y1 - 1} (608 × 608 px). Analysis was restricted to that rectangle. The two output masks retain the original 2048 × 2048 position-40 canvas and are zero outside the ROI; CSV centroids are therefore zero-based **full-frame** x/y pixel coordinates.

## Channel identification

The time-lapse acquisition metadata identified channel {movie_bf} as transmitted-light brightfield. Channel {green_channel} has a {movie_meta[green_channel]['emission_nm']:.0f}-nm emission setting and a `{movie_meta[green_channel]['lut']}` LUT, so it was used as mKikume green. Channel {red_channel} has a {movie_meta[red_channel]['emission_nm']:.0f}-nm emission setting and a `{movie_meta[red_channel]['lut']}` LUT, so it was used as photoconverted mKikume red. This spectral assignment, rather than relative brightness, determines red versus green.

In the fixed acquisition, channel {tunel_channel} has {fixed_meta[tunel_channel]['emission_nm']:.0f}-nm emission and was used as CF640R/TUNEL; channel {dapi_channel} has {fixed_meta[dapi_channel]['emission_nm']:.0f}-nm emission and was used as DAPI. Saved display LUT names do not alter those spectral assignments.

## Parasite detection and exclusions

Fluorescence was flat-background suppressed independently in red and green with a difference of Gaussians (σ=1 and 12 px), divided by a robust channel noise estimate (MAD; green={noise_g:.3f} intensity units, red={noise_r:.3f}), and combined by the larger of the two SNRs. Candidate body pixels required SNR >5 in either channel and at least 15 pixels. Large connected candidates were split only when fluorescence maxima were at least 12 px apart; otherwise they were rejected as unresolved.

The final mask contains **{int(parasite_mask.max())} clearly resolved single parasites**. **{len(clumps)} unresolved clump** was excluded: {clump_details}. This fluorescence aggregate exceeded the conservative 200-pixel single-body threshold and contained maxima too close for a defensible split; it represents at least two parasites, but their exact count was not inferred. {boundary_text} Tiny isolated SNR components below 15 pixels were treated as noise.

## Local background and proliferation index

For each parasite, `raw_red` and `raw_green` are means over exactly that parasite's labeled pixels. Local background pixels came from a 6–15 px (1.9–4.8 µm) annulus around that parasite. Pixels within 3 px of any ≥3-SNR fluorescence were removed from the annulus, as were all retained parasite pixels. The remaining local pixels were averaged separately in the raw red and raw green planes to give `background_red` and `background_green`. Corrected means and the requested index were then calculated exactly as:

`corrected red = raw red − local mean red background`

`corrected green = raw green − local mean green background`

`proliferation index = 10 − corrected red / corrected green`

## Registration

Registration was recovered automatically from the complete brightfield planes, without using `TUNEL_Overview.ROI`. SIFT descriptors were ratio-filtered, then a RANSAC similarity transform mapped fixed-image coordinates to movie coordinates:

```
x_movie = {transform[0,0]:.8f} x_fixed + {transform[0,1]:.8f} y_fixed + {transform[0,2]:.5f}
y_movie = {transform[1,0]:.8f} x_fixed + {transform[1,1]:.8f} y_fixed + {transform[1,2]:.5f}
```

This is scale={scale:.7f}, rotation={angle:.4f}°, translation=({transform[0,2]:.3f}, {transform[1,2]:.3f}) px in the above origin convention. RANSAC retained {inliers}/{total_matches} descriptor matches. Inlier landmark residuals were median {np.median(residuals):.3f} px, RMSE {np.sqrt(np.mean(residuals ** 2)):.3f} px, 95th percentile {np.percentile(residuals, 95):.3f} px, maximum {residuals.max():.3f} px. `output/_diagnostics/registration_edges.png` is an automatic visual check (movie edges red, registered fixed edges green).

## Hosts and TUNEL classification

Macrophage outlines were segmented on the last-frame brightfield with GPU Cellpose cyto3 (55-px nominal diameter, inverted brightfield). Candidate masks were retained when they contained detected infection signal, then matched one-to-one to registered DAPI nuclei; the remaining DAPI-bearing host was added by nearest unique Cellpose mask. This uses live-frame boundaries for parasite containment while transferring only nuclear TUNEL status through the registration.

DAPI nuclei were segmented automatically in the complete fixed position-40 field ({total_nuclei} nuclei). For every nucleus, mean CF640R/TUNEL was locally background-corrected with a 10–22 px nuclear annulus. A two-component Gaussian mixture (component means {tunel_means[0]:.2f} and {tunel_means[1]:.2f}; equal-posterior boundary approximately {tunel_threshold:.2f} intensity units) classified nuclei as TUNEL− or TUNEL+. No group threshold was drawn manually.

| host_cell_id | corrected nuclear TUNEL | group | retained parasites |
|---:|---:|:---:|---:|
{chr(10).join(host_lines)}

Host assignment flags: {uncertain_text}

## Group summary

{chr(10).join(median_lines)}

The requested scatter plot shows every retained parasite and a black horizontal median bar for each group.
"""
    (OUT / "report.md").write_text(report)

    print(f"Wrote {len(dataframe)} parasites, {int(host_mask_roi.max())} hosts")
    print(dataframe.to_string(index=False))


if __name__ == "__main__":
    main()
