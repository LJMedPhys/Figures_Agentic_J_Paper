#!/usr/bin/env python3
"""Minimal ImageJ .roi reader — bounds and centre, no Java.

The lab's ground truth for the auto-mode runs is positional: which parasites and
nuclei exist, and where. The workbooks carry the measurements but not the
coordinates; the `.roi` files carry the coordinates. So for a run that segmented
its own objects, the only way to ask "did it find the lab's parasite?" is to read
the ROI files.

ImageJ's format (ij.io.RoiEncoder) is a short big-endian header: magic `Iout`,
version, type, then top/left/bottom/right as int16. That is all that is needed
here — the bounding-box centre locates an object to within a pixel or two, which
is far finer than the distance at which two different parasites could be
confused.

    read(path)     -> {name, type, left, top, right, bottom, cx, cy, w, h}
    read_dir(path) -> [that, ...] for every .roi below it
"""

import os
import struct

TYPES = {0: "polygon", 1: "rect", 2: "oval", 3: "line", 4: "freeline",
         5: "polyline", 6: "noRoi", 7: "freehand", 8: "traced", 9: "angle",
         10: "point"}


def read(path):
    with open(path, "rb") as fh:
        b = fh.read(64)
    if len(b) < 18 or b[:4] != b"Iout":
        raise ValueError(f"not an ImageJ ROI: {path}")
    version, roi_type = struct.unpack(">Hb", b[4:7])
    top, left, bottom, right = struct.unpack(">hhhh", b[8:16])
    return dict(name=os.path.basename(path)[:-4] if path.endswith(".roi")
                else os.path.basename(path),
                path=path, version=version,
                type=TYPES.get(roi_type, str(roi_type)),
                left=left, top=top, right=right, bottom=bottom,
                w=right - left, h=bottom - top,
                cx=(left + right) / 2.0, cy=(top + bottom) / 2.0)


def read_dir(root, recursive=True):
    out = []
    if recursive:
        walk = ((d, fs) for d, _sub, fs in os.walk(root))
    else:
        walk = [(root, os.listdir(root))]
    for d, files in walk:
        for f in sorted(files):
            if f.lower().endswith(".roi"):
                try:
                    r = read(os.path.join(d, f))
                except (OSError, ValueError):
                    continue
                r["folder"] = os.path.relpath(d, root)
                out.append(r)
    return out


if __name__ == "__main__":
    import sys
    for arg in sys.argv[1:]:
        rois = read_dir(arg) if os.path.isdir(arg) else [read(arg)]
        print(f"{arg}: {len(rois)} ROI(s)")
        for r in rois:
            print(f"   {r['folder'] + '/' if r.get('folder', '.') != '.' else '':<12}"
                  f"{r['name']:24} {r['type']:9} "
                  f"bounds {r['left']:5d},{r['top']:5d} {r['w']:4d}x{r['h']:<4d} "
                  f"centre {r['cx']:8.1f},{r['cy']:8.1f}")
