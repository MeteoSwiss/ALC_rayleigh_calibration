# Examples

Jupyter notebooks showing how to use the `calibration` package.

Install the package first (from the repo root), then launch Jupyter from the repo root
so the notebooks find `options.json`:

```bash
uv pip install -e ".[plotting]"   # (or: pip install -e ".[plotting]")
uv pip install jupyter
jupyter lab examples/             # or: jupyter notebook
```

| Notebook | What it shows | Needs field data? |
|---|---|:--:|
| `01_calibration_quickstart.ipynb` | End-to-end Rayleigh calibration of one night: load options, define the instrument, run `calibrate_rayleigh`, read the result; plus CLI / batch usage | yes (your L1/L2 + CAMS) |
| `02_molecular_model_offline.ipynb` | Builds the molecular reference profile from the shipped US Standard Atmosphere and lists the molecular-window methods — **runs with no external data** | no |
| `03_sample_data_calibration.ipynb` | Runs Rayleigh **and** cloud calibration on the small **bundled** fixtures in `data/` (all five instrument types). CHM15k/Mini-MPL run offline; the 910 nm CL61/CL31/CL51 cells run if CAMS is present | no (uses bundled `data/`) |

The `data/` folder holds small trimmed real L1/L2 fixtures (one or two per instrument
type) — see [`data/README.md`](data/README.md). They back `03_…` and
`tests/test_sample_data.py`.

Start with `02` to confirm the install, then `03` to see real calibrations on the
bundled data; use `01` as the template once your own data paths are set.
