#!/usr/bin/env bash
# Regenerate Overview_nature.svg from ../Overview.svg and render a ~480-dpi PNG + vector PDF.
# Chromium comes from the presenton image; the host's Liberation Sans (Arial
# metrics) is mounted in because the image ships no Arial-compatible font.
set -euo pipefail
cd "$(dirname "$0")"
python3 fix_overview.py

python3 - <<'PY'
import io
svg = io.open("Overview_nature.svg", encoding="utf-8").read().split("?>", 1)[1]
io.open("_proof.html", "w", encoding="utf-8").write(
    "<!doctype html><html><head><meta charset='utf-8'><style>"
    "@page{size:180mm 121.5mm;margin:0}html,body{margin:0;padding:0}"
    "svg{display:block}</style></head><body>%s</body></html>" % svg)
PY

DOCK=(docker run --rm -v "$PWD:/w" -w /w
      -v /usr/share/fonts/truetype/liberation:/usr/share/fonts/truetype/liberation:ro
      --entrypoint chromium ghcr.io/presenton/presenton:latest
      --headless --no-sandbox --disable-gpu --hide-scrollbars --allow-file-access-from-files)
# 180 mm at 96 css-px/in = 680 px; scale 5 -> 3400 px wide (~480 dpi)
"${DOCK[@]}" --force-device-scale-factor=5 --window-size=681,460 \
      --screenshot=/w/Overview_nature.png /w/_proof.html
"${DOCK[@]}" --no-pdf-header-footer --print-to-pdf=/w/Overview_nature.pdf /w/_proof.html
rm -f _proof.html
echo "-> Overview_nature.svg / .png / .pdf"
