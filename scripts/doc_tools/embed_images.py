"""Embed local image references in a markdown report as base64 data URIs, producing a
self-contained <name>_embedded.md (images always render — portable to any viewer / Word / PDF).

Usage: python embed_images.py report1.md [report2.md ...]
"""
import re, base64, sys
from pathlib import Path

_MIME = {"png": "png", "jpg": "jpeg", "jpeg": "jpeg", "gif": "gif", "svg": "svg+xml"}


def embed(md_path):
    md = Path(md_path)
    text = md.read_text(encoding="utf-8")
    base = md.parent
    n_ok = n_miss = 0

    def repl(m):
        nonlocal n_ok, n_miss
        alt, path = m.group(1), m.group(2).strip()
        if path.startswith(("data:", "http://", "https://")):
            return m.group(0)
        p = (base / path)
        if not p.exists():
            n_miss += 1
            return m.group(0)
        mime = _MIME.get(p.suffix.lstrip(".").lower(), "png")
        b = base64.b64encode(p.read_bytes()).decode("ascii")
        n_ok += 1
        return f"![{alt}](data:image/{mime};base64,{b})"

    out = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", repl, text)
    dest = md.with_name(md.stem + "_embedded.md")
    dest.write_text(out, encoding="utf-8")
    print(f"  {md.name}: embedded {n_ok} images ({n_miss} missing) -> {dest.name} ({len(out)/1e6:.1f} MB)")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        embed(a)
