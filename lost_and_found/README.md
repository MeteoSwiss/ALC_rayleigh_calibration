# lost_and_found/

One-off R&D scripts, scratch / debug helpers and old run logs, kept for reference.

These are **not maintained** and are not part of the supported package or workflow. They
may hardcode paths, import sibling scripts, or depend on transient data. If something
here turns out to be genuinely reusable, promote it to `scripts/` or `validation/`.

Contents:

- **Ad-hoc sensitivity probes** — `run_cl61_rayleigh_params.py`, `run_cl61_rayleigh_nightwindow.py`,
  `run_cl61_sza_test.py`, `sweep_gates_payerne.py`, `run_camsmol_payerne.py`, `run_rayleigh_diag_dates.py`
- **Night-selection helpers** — `select_nights.py`, `select_nights_site.py`
- **Quick plots / counts** — `plot_cl_timeseries.py`, `count_flags_payerne.py`
- **Profiling / inspection scratch** — `_prof.py`, `_res.py`, `_make_lindenberg_report_figs.py`
  (the last one reuses helpers from `scripts/run_lindenberg_cl61_cal.py`)
- **Old run logs** — `*.log` (git-ignored)
