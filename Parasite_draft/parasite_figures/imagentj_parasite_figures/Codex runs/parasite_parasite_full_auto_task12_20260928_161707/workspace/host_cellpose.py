#!/usr/bin/env python3
"""Headless Cellpose helper used by analyze_parasite_task.py."""

import argparse
import numpy as np
import tifffile
from cellpose import core, models


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tif")
    parser.add_argument("output_npy")
    args = parser.parse_args()

    image = tifffile.imread(args.input_tif)
    model = models.CellposeModel(gpu=core.use_gpu(), model_type="cyto3")
    normalize = {
        "lowhigh": None,
        "percentile": [1, 99],
        "normalize": True,
        "norm3D": False,
        "sharpen_radius": 0,
        "smooth_radius": 0,
        "tile_norm_blocksize": 0,
        "tile_norm_smooth3D": 1,
        "invert": False,
    }
    masks, _, _ = model.eval(
        image,
        diameter=90,
        channels=[0, 0],
        normalize=normalize,
        flow_threshold=0.4,
        cellprob_threshold=-1.0,
    )
    np.save(args.output_npy, masks.astype(np.uint16))


if __name__ == "__main__":
    main()
