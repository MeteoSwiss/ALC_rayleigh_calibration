"""
embed_report.py — produce a self-contained <name>_embedded.md from a report whose images are relative
links into the external figure tree, by inlining each PNG as a base64 data URI.

Usage:  python validation/embed_report.py <report.md> [figs_base_dir]
        figs_base_dir defaults to the project figure root (where 'figs_paper_validation/...' resolves).
"""
from __future__ import annotations
import base64
import re
import sys
from pathlib import Path

DEFAULT_BASE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
IMG = re.compile(r"!\[([^\]]*)\]\(([^)]+\.png)\)")


def main():
    report = Path(sys.argv[1])
    base = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_BASE
    text = report.read_text(encoding="utf-8")
    missing = []

    def repl(m):
        caption, rel = m.group(1), m.group(2)
        if rel.startswith("data:"):
            return m.group(0)
        png = (base / rel)
        if not png.is_file():
            png = report.parent / rel
        if not png.is_file():
            missing.append(rel)
            return m.group(0)
        b64 = base64.b64encode(png.read_bytes()).decode("ascii")
        return f"![{caption}](data:image/png;base64,{b64})"

    out = IMG.sub(repl, text)
    dst = report                      # embed in place -> a single report file with figures
    dst.write_text(out, encoding="utf-8")
    print(f"embedded figures in place: {dst}  ({dst.stat().st_size/1e6:.2f} MB)")
    if missing:
        print("  WARNING missing figures:", *missing, sep="\n   ")


if __name__ == "__main__":
    main()
