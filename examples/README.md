# Examples

Jupyter notebooks showing how to use the `rayleigh_calibration` package.

Install the package first (from the repo root), then launch Jupyter from the repo root
so the notebooks find `options.json`:

```bash
pip install -e ".[plotting]"
pip install jupyter
jupyter lab examples/        # or: jupyter notebook
```

| Notebook | What it shows | Needs field data? |
|---|---|:--:|
| `01_rayleigh_calibration_quickstart.ipynb` | End-to-end Rayleigh calibration of one night: load options, define the instrument, run `calibrate_rayleigh`, read the result; plus CLI / batch usage | yes (L1/L2 + CAMS) |
| `02_molecular_model_offline.ipynb` | Builds the molecular reference profile from the shipped US Standard Atmosphere and lists the molecular-window methods — **runs with no external data** | no |

Start with `02` to confirm the install works, then `01` once your data paths are set.
