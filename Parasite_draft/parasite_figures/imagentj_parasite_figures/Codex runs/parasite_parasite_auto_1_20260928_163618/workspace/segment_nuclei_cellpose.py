#!/usr/bin/env python3
"""Cellpose-SAM helper, run with the dedicated cellpose4 environment."""

from __future__ import annotations

import argparse

import numpy as np
import tifffile
from cellpose import models


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tif")
    parser.add_argument("output_tif")
    args = parser.parse_args()

    image = tifffile.imread(args.input_tif)
    model = models.CellposeModel(gpu=True)
    masks, _, _ = model.eval(image, diameter=35)
    tifffile.imwrite(args.output_tif, np.asarray(masks, dtype=np.uint16))


if __name__ == "__main__":
    main()
