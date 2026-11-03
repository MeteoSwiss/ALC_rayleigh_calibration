"""Slide 8 of the Annual Meeting 2026 deck, v2.3 edition ("Which method for which instrument"): the slide
redraw of paper Fig. 11 (presentation/slide8_methods.py) fed with the manuscript v1.0 data, i.e. the forced
cloud + Rayleigh run of the 33 test instruments under the calibration v2.3 release configuration
(validation/paper/fig_v10_method_comparison.py: eprof_v2.3, uncertainty budget, polar fallback, eta tables v2,
February - June 2026, C:/DATA/Projects/202611_v23_release/method_test_v23).

The slide module is loaded from its file and its paper module swapped for the v1.0 one, so the panels, fonts
and statistics are exactly those of the 2 November deck.

Run from the repo root:  python presentation/slide8_methods_v23.py
Output: presentation/figs_annual2026/fig_methods_v23_slide.png (+ the per-type summary JSON)
"""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from validation.paper import fig_v10_method_comparison as M10  # noqa: E402

spec = importlib.util.spec_from_file_location("slide8_methods", REPO / "presentation" / "slide8_methods.py")
S8 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S8)
S8.M = M10
S8.OUT_PNG = REPO / "presentation" / "figs_annual2026" / "fig_methods_v23_slide.png"

if __name__ == "__main__":
    S8.main()
