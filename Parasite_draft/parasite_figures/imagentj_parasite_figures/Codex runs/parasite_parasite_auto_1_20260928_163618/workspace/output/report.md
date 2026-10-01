# TUNEL signal in nuclei — Position040

## Source image, channels, and ROI check

I analyzed series **Position040** (series index 39) from `20210728_TUNEL_staining.lif`. The Leica hardware metadata identifies channel 0 as transmitted-light brightfield, channel 1 as the 635-nm-excited/695-nm-emission acquisition (CF640R; TUNEL), and channel 2 as the 385-nm-excited/425-nm-emission acquisition (DAPI). Thus channel 1 was measured as TUNEL and channel 2 was used for nuclei.

Both fluorescence channels are 2048 × 2048 pixels and match one another exactly. The ImageJ ROI is a rectangle x=565…1172, y=952…1559 (width 608, height 608 pixels; zero-based coordinates). An ImageJ `.roi` stores its bounding coordinates, not the dimensions of its source image, so literal equality between ROI size and image size is not expected: this ROI deliberately specifies a smaller analysis region. Its complete bounding box is valid within the 2048 × 2048 source, and the extracted region in both channels is 608 × 608 pixels. No resize, registration, or coordinate correction was needed. Pixel size from the LIF metadata is 0.3223 µm/pixel. The delivered label mask is 2048 × 2048 so that it aligns directly with Position040; pixels outside the supplied ROI are zero.

## Nucleus segmentation

Segmentation was fully automatic and used only the DAPI image inside the supplied ROI. I used the pretrained **Cellpose 4 Cellpose-SAM (`cpsam`)** model with GPU inference and a 35-pixel target diameter (11.3 µm at this image's sampling). This model resolved the one touching/overlapping nuclear pair in the ROI. I then retained plausible nuclear objects of 90–2500 px and excluded partial objects touching the analysis-ROI boundary. IDs were reassigned consecutively from top to bottom. This produced **7 nuclei**. Centroids in the CSV are x/y coordinates in the full Position040 image, not ROI-local coordinates.

## TUNEL measurement and local background

`raw_tunel_intensity` is the mean channel-1 intensity over every pixel of a nucleus. For each nucleus independently, I formed an automatically generated local annulus extending 6–18 px outside that nucleus. All segmented nuclear pixels were excluded, including neighboring nuclei. To keep isolated bright extracellular structures or hot pixels from biasing the estimate, I discarded the brightest 20% of the remaining annular pixels and averaged the lower 80%. This is the reported `background_intensity`. Therefore every background was sampled immediately next to the measured nucleus, without hand-drawn regions. `corrected_tunel_intensity = raw_tunel_intensity - background_intensity`.

## TUNEL-positive cut

The cut was selected objectively from this ROI rather than supplied in advance: **two-class Otsu on log1p-transformed corrected intensities, with the final cut centered in log space between the brightest low-class nucleus and the dimmest high-class nucleus**, giving **16.603 camera units**. I used the log transform because the bright-positive population is right-skewed, and centered the final cut in the observed empty interval so that it does not sit artificially close to either measured nucleus. Values were classified TUNEL+ only when strictly above the cut. This yielded **3 TUNEL+** and **4 TUNEL-** nuclei. Median corrected intensities were 13.100 (TUNEL-) and 21.628 (TUNEL+).

## Resolution and caveats

All requested quantities were resolved. The classification is an image-intensity phenotype, not an independent biological validation of apoptosis. Nuclei cut by the supplied ROI boundary were deliberately excluded because neither their area nor mean intensity could be measured completely.
