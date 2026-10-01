# Position 40, last-frame parasite proliferation analysis

## Scope and coordinate convention

The supplied `LastFrame_Overview.ROI` was read directly as x=537:1144, y=849:1456 (608 × 608 px). Analysis was restricted to that rectangle. The two output masks retain the original 2048 × 2048 position-40 canvas and are zero outside the ROI; CSV centroids are therefore zero-based **full-frame** x/y pixel coordinates.

## Channel identification

The time-lapse acquisition metadata identified channel 0 as transmitted-light brightfield. Channel 1 has a 525-nm emission setting and a `Green` LUT, so it was used as mKikume green. Channel 2 has a 600-nm emission setting and a `Red` LUT, so it was used as photoconverted mKikume red. This spectral assignment, rather than relative brightness, determines red versus green.

In the fixed acquisition, channel 1 has 695-nm emission and was used as CF640R/TUNEL; channel 2 has 425-nm emission and was used as DAPI. Saved display LUT names do not alter those spectral assignments.

## Parasite detection and exclusions

Fluorescence was flat-background suppressed independently in red and green with a difference of Gaussians (σ=1 and 12 px), divided by a robust channel noise estimate (MAD; green=5.303 intensity units, red=2.235), and combined by the larger of the two SNRs. Candidate body pixels required SNR >5 in either channel and at least 15 pixels. Large connected candidates were split only when fluorescence maxima were at least 12 px apart; otherwise they were rejected as unresolved.

The final mask contains **11 clearly resolved single parasites**. **1 unresolved clump** was excluded: bbox (x=227:250, y=299:318), area=304 px, peaks=2, minimum peak spacing=9.2 px. This fluorescence aggregate exceeded the conservative 200-pixel single-body threshold and contained maxima too close for a defensible split; it represents at least two parasites, but their exact count was not inferred. In addition, 1 strong fluorescence object at the left ROI boundary was rejected as boundary-truncated (not counted as a clump), because its full shape could not be evaluated. Tiny isolated SNR components below 15 pixels were treated as noise.

## Local background and proliferation index

For each parasite, `raw_red` and `raw_green` are means over exactly that parasite's labeled pixels. Local background pixels came from a 6–15 px (1.9–4.8 µm) annulus around that parasite. Pixels within 3 px of any ≥3-SNR fluorescence were removed from the annulus, as were all retained parasite pixels. The remaining local pixels were averaged separately in the raw red and raw green planes to give `background_red` and `background_green`. Corrected means and the requested index were then calculated exactly as:

`corrected red = raw red − local mean red background`

`corrected green = raw green − local mean green background`

`proliferation index = 10 − corrected red / corrected green`

## Registration

Registration was recovered automatically from the complete brightfield planes, without using `TUNEL_Overview.ROI`. SIFT descriptors were ratio-filtered, then a RANSAC similarity transform mapped fixed-image coordinates to movie coordinates:

```
x_movie = 0.99904857 x_fixed + -0.02989639 y_fixed + 31.08984
y_movie = 0.02989639 x_fixed + 0.99904857 y_fixed + -126.10774
```

This is scale=0.9994958, rotation=1.7141°, translation=(31.090, -126.108) px in the above origin convention. RANSAC retained 367/372 descriptor matches. Inlier landmark residuals were median 0.627 px, RMSE 0.772 px, 95th percentile 1.358 px, maximum 1.977 px. `output/_diagnostics/registration_edges.png` is an automatic visual check (movie edges red, registered fixed edges green).

## Hosts and TUNEL classification

Macrophage outlines were segmented on the last-frame brightfield with GPU Cellpose cyto3 (55-px nominal diameter, inverted brightfield). Candidate masks were retained when they contained detected infection signal, then matched one-to-one to registered DAPI nuclei; the remaining DAPI-bearing host was added by nearest unique Cellpose mask. This uses live-frame boundaries for parasite containment while transferring only nuclear TUNEL status through the registration.

DAPI nuclei were segmented automatically in the complete fixed position-40 field (51 nuclei). For every nucleus, mean CF640R/TUNEL was locally background-corrected with a 10–22 px nuclear annulus. A two-component Gaussian mixture (component means 5.37 and 19.53; equal-posterior boundary approximately 13.66 intensity units) classified nuclei as TUNEL− or TUNEL+. No group threshold was drawn manually.

| host_cell_id | corrected nuclear TUNEL | group | retained parasites |
|---:|---:|:---:|---:|
| 1 | 3.85 | TUNEL- | 1 |
| 2 | 5.83 | TUNEL- | 0 |
| 3 | 15.30 | TUNEL+ | 1 |
| 4 | 2.50 | TUNEL- | 3 |
| 5 | 4.78 | TUNEL- | 3 |
| 6 | 11.07 | TUNEL- | 3 |

Host assignment flags: None; every retained parasite overlapped one retained host mask by at least 50% of its pixels.

## Group summary

- TUNEL-: n=10, median=7.340
- TUNEL+: n=1, median=9.659

The requested scatter plot shows every retained parasite and a black horizontal median bar for each group.
