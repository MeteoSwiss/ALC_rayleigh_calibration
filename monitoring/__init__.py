"""Daily monitoring dashboard for the E-PROFILE Rayleigh calibration outputs.

A lightweight, static-site generator: it indexes the per-station calibration CSVs into
a single SQLite file and renders an HTML summary plus one page per instrument. No server
or database process is required -- the output is a folder of files served from any web dir.

Pipeline:  build_index() -> SQLite  ->  build_site() -> HTML

See scripts/build_dashboard.py for the CLI entry point.
"""
