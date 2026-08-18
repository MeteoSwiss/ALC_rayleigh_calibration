# -*- coding: utf-8 -*-
"""Assemble un npz de dark ESTIME par station, consommable par ALC_DARK_PROFILE.

Pour chaque ident demande, prend la reconstruction template du scan nuits-claires
(theta * forme capot du type, cle `<ident>_b_tpl` du npz par flux) -- PAS le profil libre,
prouve contamine en amplitude. La cle de sortie est `<ident>_b_rcs`, celle que lit
`calibration.rayleigh.calibration._dark_on_grid`.

Run : python rayleigh_availability/build_dark_npz.py --wmo 0-20000-0-06240 --idents A,B,C,D
      python rayleigh_availability/build_dark_npz.py --wmo 0-20000-0-10393 --idents 0,C
Sortie : <DATA>/dark_clearsky/dark_est_<wmo>.npz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/dark_clearsky")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wmo", required=True)
    ap.add_argument("--idents", required=True, help="liste separee par des virgules")
    args = ap.parse_args()

    out = {}
    for ident in args.idents.split(","):
        ident = ident.strip()
        src = DATA / f"{args.wmo}_{ident}.npz"
        if not src.exists():
            print(f"  {ident}: npz absent ({src.name}), saute")
            continue
        z = np.load(src)
        b = np.asarray(z[f"{ident}_b_tpl"], float)
        rng = np.asarray(z[f"{ident}_range"], float)
        theta = float(z["theta"]) if "theta" in z else np.nan
        if not np.any(np.isfinite(b)) or not np.isfinite(theta):
            print(f"  {ident}: pas de reconstruction template (theta={theta}), saute")
            continue
        out[f"{ident}_range"] = rng
        out[f"{ident}_b_rcs"] = np.nan_to_num(b)     # la cle lue par _dark_on_grid
        print(f"  {ident}: theta={theta:+.3f}, b range {np.nanmin(b):+.3e}..{np.nanmax(b):+.3e}")
    if not out:
        print("rien a ecrire")
        return
    dst = DATA / f"dark_est_{args.wmo}.npz"
    np.savez(dst, **out)
    print(f"-> {dst}")


if __name__ == "__main__":
    main()
