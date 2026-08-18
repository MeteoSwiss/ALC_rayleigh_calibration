#!/usr/bin/env python3
"""Thin shim: the daily-calibration panel is PRODUCTION code and lives in monitoring/panel.py.

It started life as this mockup; the production renderer then imported production behaviour from a
file named mockup_*, which is how two generations of page code ended up sharing one site (see
doc/reports/12_dashboard_professionalization_plan.md). This shim keeps old command lines working.
"""
from monitoring.panel import *            # noqa: F401,F403
from monitoring.panel import main

if __name__ == "__main__":
    main()
