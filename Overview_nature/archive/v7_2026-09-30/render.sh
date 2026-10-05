#!/usr/bin/env bash
# Regenerate Overview_nature.svg (panel a alone) and Fig1_overview_modes.svg (a + b mode
# tree) and render each to a ~480-dpi PNG + vector PDF at its true mm size.
# Chromium comes from the presenton image; the host's Liberation Sans (Arial
# metrics) is mounted in because the image ships no Arial-compatible font.
# Archive the previous outputs into archive/vN_YYYY-MM-DD/ before re-running.
set -euo pipefail
cd "$(dirname "$0")"
python3 fix_overview.py
python3 make_fig1.py

DOCK=(docker run --rm -v "$PWD:/w" -w /w
      -v /usr/share/fonts/truetype/liberation:/usr/share/fonts/truetype/liberation:ro
      --entrypoint chromium ghcr.io/presenton/presenton:latest
      --headless --no-sandbox --disable-gpu --hide-scrollbars --allow-file-access-from-files)

render() {
  local name=$1 W H pw ph
  read -r W H < <(python3 - "$name.svg" <<'PY'
import re, sys
head = open(sys.argv[1], encoding="utf-8").read(4000)
tag = re.search(r"<svg\b[^>]*>", head, re.S).group(0)
print(re.search(r'\bwidth="([0-9.]+)mm"', tag).group(1), re.search(r'\bheight="([0-9.]+)mm"', tag).group(1))
PY
)
  python3 - "$name" "$W" "$H" <<'PY'
import io, sys
name, w, h = sys.argv[1:4]
svg = io.open(name + ".svg", encoding="utf-8").read().split("?>", 1)[1]
io.open("_proof.html", "w", encoding="utf-8").write(
    "<!doctype html><html><head><meta charset='utf-8'><style>"
    f"@page{{size:{w}mm {h}mm;margin:0}}html,body{{margin:0;padding:0}}"
    "svg{display:block}</style></head><body>%s</body></html>" % svg)
PY
  # 96 css-px per inch; scale 5 -> ~480 dpi
  pw=$(python3 -c "import math;print(math.ceil($W/25.4*96)+1)")
  ph=$(python3 -c "import math;print(math.ceil($H/25.4*96)+1)")
  "${DOCK[@]}" --force-device-scale-factor=5 --window-size="$pw,$ph" \
        --screenshot="/w/$name.png" /w/_proof.html 2>/dev/null
  "${DOCK[@]}" --no-pdf-header-footer --print-to-pdf="/w/$name.pdf" /w/_proof.html 2>/dev/null
  rm -f _proof.html
  echo "-> $name.svg / .png / .pdf  (${W} x ${H} mm)"
}
render Overview_nature
render Fig1_overview_modes
