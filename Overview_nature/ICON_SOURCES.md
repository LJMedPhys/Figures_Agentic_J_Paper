# Fig. 1 — sources of all icons and logos

Applies to `Fig1_overview_modes.svg` / `Overview_nature.svg` from v10 (2026-10-01) on (logos: from v13, 2026-10-02).
Everything not listed here (boxes, arrows, text, the person icons in panel b) is
drawn by the authors.

## Icons — Tabler Icons (MIT License, no attribution required)

Source: Tabler Icons by Paweł Kuna, <https://tabler.io/icons>,
<https://github.com/tabler/tabler-icons>, release **v3.48.0**, outline style.
License: MIT — copy in `icons/tabler/LICENSE`. The SVG files used are in `icons/tabler/`.
Line width reduced from 2 to 1.25 (Tabler units) and colour set to black; shapes unchanged.

| Figure element | Tabler icon | File (pinned) |
|---|---|---|
| Watchdog | `user-shield` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/user-shield.svg> |
| Visual Judge | `eye` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/eye.svg> |
| Quality Assurance | `zoom-check` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/zoom-check.svg> |
| Discovery | `puzzle` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/puzzle.svg> |
| Skills: Plotting Guidelines | `chart-dots` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/chart-dots.svg> |
| Skills: Statistics | `chart-histogram` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/chart-histogram.svg> |
| Skills: Publications Standards | `building-bank` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/building-bank.svg> |
| Recipe & Mistakes DB | `database` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/database.svg> |
| Librarian | `book` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/book.svg> |
| Tools: File Tools | `folder-open` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/folder-open.svg> |
| Tools: Inspect UI | `device-desktop-analytics` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/device-desktop-analytics.svg> |
| Tools: Check Environment | `list-check` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/list-check.svg> |
| Tools: Collect Metadata | `database-cog` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/database-cog.svg> |
| Tools: Script Execution | `code` | <https://raw.githubusercontent.com/tabler/tabler-icons/v3.48.0/icons/outline/code.svg> |

## Logos — official files, unmodified, to refer to the software

Checked 2026-10-02 against each project's own files (rendered side by side).
Files are kept in `icons/logos/`; Fiji and Python are embedded byte for byte as SVG,
at their native aspect ratio (`preserveAspectRatio="xMidYMid meet"`).

| Logo | Used in | File (official source) | Terms (as of 2026-10-02) |
|---|---|---|---|
| Fiji | Browser, Fiji Coder / Debugger, Skills (Fiji plugins) | `fiji.svg` (90 x 90, incl. its own soft shadow) from <https://raw.githubusercontent.com/imagej/imagej.github.io/main/media/icons/fiji.svg> — byte-identical to `images/fiji-logo-1.0.svg` in <https://github.com/fiji/fiji> (there since 2009; Inkscape export path `/home/albert/webs/fiji/`, i.e. presumably drawn by Albert Cardona) | **CC BY 3.0** as hosted on imagej.net: the site's source repo <https://github.com/imagej/imagej.github.io> carries `LICENSE.txt` = CC BY 3.0 Unported (since the first commit, 2020-03-09) and nothing exempts `media/`. The same file is also in the GPL-3 fiji/fiji repo. imagej.net itself has **no** copyright, logo or trademark statement (footer, /licensing, /software/fiji, /contribute/governance all checked 2026-10-02), so no ™ is required. Community practice (image.sc, <https://forum.image.sc/t/license-usage-of-the-fiji-logo/20357>): use **unmodified** to refer to Fiji, following the intent of the PSF policy. Wikimedia Commons: below threshold of originality — <https://commons.wikimedia.org/wiki/File:FIJI_(software)_Logo.svg>. |
| napari | Browser, Skills (napari plugins) | embedded PNG, 824 x 824, pixel-identical in appearance to `gradient.svg` from <https://github.com/napari/napari-resources/tree/main/src/napari_resources/resources/logos/variants> | **CC BY-NC-ND 4.0** (napari-resources README). Guidelines explicitly allow the logo "to describe factual use of the software in your own projects, for example in papers, figures, posters"; must **not** imply endorsement by napari. BY → credit required; ND → no alterations; spelled **napari** (lower case). |
| Python | Python Coder, Skills (Python packages) | `python-logo-only.svg` (the official "two snakes" device, incl. its soft shadow) from <https://www.python.org/community/logos/> → <https://s3.dualstack.us-east-2.amazonaws.com/pythondotorg-assets/media/files/python-logo-only.svg> | PSF policy <https://www.python.org/psf/trademarks/>: "The Python logos are not currently registered", so they carry **™** (not ®), which "must always be included with the logo" → a 5-pt ™ sits at the top right of both logos. Unaltered colours and shapes only. ® on the *word* "Python" is waived in academic papers. |

### Problems found and fixed in v13

- Fiji: the old embedded copy was a blurry 90 x 90 px raster of `fiji.svg` (~280 dpi at
  print size, below the 300 dpi minimum). It is now the vector original.
- Fiji (coder box, Skills) and Python (both places) were squashed horizontally by 9.1 %
  (`preserveAspectRatio="none"` on a non-square box), which counts as an alteration
  under every one of the three policies. All logos are now at their native aspect ratio.
- Python: the ™ was missing; the old copy also had no shadow (origin unknown), so it is
  replaced by the PSF's own logo-only file.
- napari: written "Napari" in the figure; now "napari".

### Open points for the authors

- Nature Methods publishes figures under the article's licence (CC BY for open access).
  The napari logo is CC BY-NC-ND, so the legend should name its licence and state that it is
  excluded from the article's licence. Before submission, consider getting written OK from
  the napari Steering Council (cheap insurance, because NC sits uneasily with a commercial publisher).
- The image.sc threads could not be re-read on 2026-10-02 (the forum blocks automated
  access); the summaries above come from the earlier version of this file and from search snippets.

## Suggested credit line for the Fig. 1 legend

> Icons: Tabler Icons (MIT licence). Python™ logo: "Python" and the Python logos are
> trademarks or registered trademarks of the Python Software Foundation. Fiji logo:
> imagej.net / Fiji contributors, CC BY 3.0. napari logo:
> © napari contributors, CC BY-NC-ND 4.0 (https://github.com/napari/napari-resources),
> not covered by this article's licence. The Fiji and napari logos are used unmodified
> to refer to the respective software and do not imply endorsement.

## History

- v1–v9: 14 icons traced (Inkscape bitmap trace) from PNGs of unknown origin.
- v10: all 14 replaced by Tabler icons (Discovery = `brain`); logos unchanged.
- v11: Discovery switched to `puzzle` at the authors' request.
- v12: subagent icons/logos aligned in one column.
- v13 (2026-10-02): Fiji + Python replaced by the official SVGs at native aspect
  (were blurry / squashed 9 %); ™ added to both Python logos; "Napari" -> "napari".
