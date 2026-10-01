#!/usr/bin/env python3
"""Cellpose stage for analyze_parasites.py (run in the cellpose environment)."""

import argparse
from pathlib import Path

import numpy as np
import tifffile
from cellpose import models
from skimage.exposure import equalize_adapthist


def robust_u8(image):
    lo, hi = np.percentile(image, [1, 99])
    scaled = np.clip((image.astype(float) - lo) / max(hi - lo, 1), 0, 1)
    return np.uint8(255 * scaled)


parser = argparse.ArgumentParser()
parser.add_argument("last_frame_roi_tif", type=Path)
parser.add_argument("output_mask_tif", type=Path)
args = parser.parse_args()

planes = tifffile.imread(args.last_frame_roi_tif)
bf = robust_u8(planes[0])
bf = np.uint8(255 * equalize_adapthist(bf, clip_limit=0.015))

model = models.CellposeModel(gpu=True, model_type="cyto3")
mask, _, _ = model.eval(
    bf,
    channels=[0, 0],
    diameter=55,
    invert=True,
    flow_threshold=0.5,
    cellprob_threshold=-0.2,
    min_size=150,
)
tifffile.imwrite(args.output_mask_tif, mask.astype(np.uint16))
