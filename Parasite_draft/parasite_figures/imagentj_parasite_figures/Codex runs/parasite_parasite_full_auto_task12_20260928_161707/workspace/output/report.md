# Position 40: TUNEL and parasite proliferation analysis

## Inputs, channels, and ROI checks

Both position-40 source series are 2048 x 2048 pixels. The stain ROI is `(x=565:1173, y=952:1560)`; the last-frame ROI is `(x=537:1145, y=849:1457)`. Both are valid 608 x 608 rectangles entirely inside their source images. ImageJ rectangle ROI files do not store the source canvas width/height, so an independent ROI-canvas dimension cannot be compared; the bounds, equal ROI sizes, and cell content were checked against the native 2048 x 2048 sources. No rescaling or ROI correction was needed.

The stained LIF acquisition settings were 695 nm for channel 1 and 425 nm for channel 2. Channel 1 is therefore the far-red CF640R/TUNEL acquisition and channel 2 is DAPI. The saved viewer LUT names `['Gray', 'Blue', 'Cyan']` are misleading relative to the stain identities; physical acquisition settings and image morphology were used. In the movie, the LIF LUT metadata identifies channel 1 as green mKikume and channel 2 as red mKikume (`['Gray', 'Green', 'Red']`).

## Nuclear TUNEL analysis

DAPI was Gaussian-smoothed (sigma 1.2 px), thresholded at ROI median + 8 robust MAD (threshold 196.16), closed by 2 px, and filtered to 80-2500 px. This yielded 6 nuclei. IDs and CSV centroids are in full stained-image coordinates; the TIFF is a full 2048 x 2048 label image with zero outside the ROI.

For each DAPI nucleus, the TUNEL value is the mean of its pixels. Local background was measured in a 6-22 px annulus around that nucleus, excluding all nuclei plus a 2 px guard. Pixels outside median +/- 3 robust MAD within the annulus were rejected to avoid bright debris, and the mean of the remaining local pixels was subtracted. Thus every nucleus has its own nearby background estimate.

The corrected nuclear values were 4.93, 5.05, 6.01, 7.68, 10.79, 13.21 a.u. The largest gap was between 7.68 and 10.79; I placed the cut at their midpoint, **9.23 a.u.** Values above it are TUNEL+, values below it TUNEL-. This gives 2 TUNEL+ and 4 TUNEL- nuclei. This is an image-internal intensity cut; with only six nuclei it should not be treated as a universal assay threshold.

## Registration and macrophage hosts

Registration used only brightfield: SIFT feature matches and RANSAC supplied the initialization, followed by rigid ECC refinement on background-subtracted images. The local stain-ROI to movie-ROI affine matrix `[x'; y'] = A[x; y; 1]` is `[0.999505, -0.031470, 30.137085; 0.031470, 0.999505, -6.334007]`; the equivalent full-image matrix is `[0.999505, -0.031470, 32.376676; 0.031470, 0.999505, -126.643195]`. It is a 1.803-degree rotation plus translation, with no fitted scale or shear. There were 9 RANSAC inliers among 76 ratio-tested matches; their post-refinement median residual was **1.63 px** (mean 2.00 px), and ECC correlation was 0.263. Registered cell and nuclear landmarks overlay correctly.

Macrophages were segmented headlessly from last-frame brightfield with pretrained Cellpose `cyto3` (diameter 90 px, cell-probability threshold -1), then restricted to 2500-15000 px objects. It found 6 hosts. Each registered DAPI centroid matched a distinct host; five centroids were inside their host mask and one required a 6.63 px nearest-boundary association (about 2.12 micrometres), consistent with the thin nuclear edge omitted by the brightfield cell mask. Host IDs in `host_mask.tif` are the IDs used in `parasite_index.csv`.

## Parasites, clumps, and intensity calculation

The movie metadata identifies channel 2 as red and channel 1 as green. Parasites were detected in red because its local signal-to-background was stronger. I subtracted a sigma-12-px Gaussian background from a sigma-1-px image, thresholded at median + 6 robust MAD, and used persistent peaks above median + 10 robust MAD as watershed seeds. Contacts with two spatially separate persistent peaks were split. A candidate with area greater than 1.8 times the single-object median and no separable second peak was classified as an unresolved clump.

I retained 10 clearly resolved parasites and excluded 1 unresolved clump candidate(s): area 284 px at local centroid (x=237.8, y=308.8). Exclusion was by size plus lack of a resolvable internal intensity valley/second peak, not by host group.

For each parasite and each color separately, background was the sigma-clipped mean in a nearby 4-15 px annulus, excluding all detected parasite signal plus a 2 px guard. The reported raw values are within-mask means; corrected values are raw minus the corresponding local background. The proliferation index was calculated exactly as `10 - corrected_red / corrected_green`.

All 10 retained parasites fell inside a Cellpose host mask and were assigned confidently; none was guessed or left unassigned. Parasite group is the registered nuclear TUNEL class of that host. Parasite ID(s) 7 had corrected green below 3 a.u.; the requested ratio is reported, but it is denominator-sensitive and should be interpreted cautiously.

## Consistency checks and unresolved points

`tunel_nuclei_mask.tif` has 6 nonzero IDs and `tunel_nuclei.csv` has 6 rows. `parasite_mask.tif` has 10 nonzero IDs and `parasite_index.csv` has 10 rows. No manual ROI, click, or interactive correction was used. The only material uncertainty is the near-background green denominator noted above and the biological generalizability of a six-nucleus TUNEL cut.
