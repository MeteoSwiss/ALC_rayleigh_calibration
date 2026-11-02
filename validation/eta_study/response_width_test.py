"""Diagnostic for step 2: which effective instrument response width do the observed in-cloud
shapes want, and how does the shape-derived alpha move with it?

The step-2 fit (fit_alpha_shape.py) prefers the widest tabulated pulse (160 ns) even for the
110-ns CL31/CL51, i.e. the observed peak is broader than PVC + nominal Gaussian pulse + gate
integration. Here the homogeneous PVC profiles (forward_lut_ext.npz, beta_1m, a_G 3.3 um) are
convolved with Gaussian responses of 110..400 ns FWHM, gate-sampled exactly as in the fit (5 phases,
peak-aligned, unit-sum on peak-60..peak+300 m), and fitted to a random sample of valid profiles per
type. Also reports the tail-slope estimator alpha_k = k_app / (dk/dalpha of the model), which is
insensitive to the response width (convolution preserves an exponential slope).

Run:  python validation/eta_study/response_width_test.py
Output: C:/DATA/Projects/202611_eta_study/step2/response_width_test.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_alpha_shape as F  # noqa: E402

RESP_NS = (110.0, 160.0, 200.0, 250.0, 300.0, 400.0)
N_SAMPLE = 20000
C = 299_792_458.0


def conv_gauss(b, ns):
    sig = C * ns * 1e-9 / 2.0 / 2.354820045
    k = np.arange(-int(5 * sig), int(5 * sig) + 1)
    w = np.exp(-0.5 * (k / sig) ** 2)
    w /= w.sum()
    out = np.apply_along_axis(lambda x: np.convolve(x, w, mode="same"), -1, b)
    return out


def main():
    import pandas as pd
    L = F.lut()
    z = L["z_rel_m"]
    ag_i = np.searchsorted(L["aG_um"], 3.3)                        # between 3.0 and 3.5
    wa = (3.3 - L["aG_um"][ag_i - 1]) / (L["aG_um"][ag_i] - L["aG_um"][ag_i - 1])
    df = pd.concat([pd.read_parquet(p, columns=["key", "itype", "grid", "valid", "cbh_m", "time_ns", "k_app"])
                    for p in sorted((F.OUT / "profiles").glob("*_working.parquet"))], ignore_index=True)
    rng = np.random.default_rng(0)
    rows = []
    for it in ("CL31", "CL51", "CL61"):
        ti = list(L["types"]).index(it)
        b1 = L["beta_1m"][ti]                                       # (alpha16, aG8, cbh15, z)
        b33 = b1[:, ag_i - 1] * (1 - wa) + b1[:, ag_i] * wa         # (alpha, cbh, z)
        dg = float(L[f"{it}_work_gate_m"])
        ks = F.gate_offsets(dg)
        zg = ks * dg
        # sample of valid profiles at the type's standard gate
        sub = df[(df.itype == it) & df.valid]
        keys = []
        for k in sub.key.unique():
            zz = np.load(F.SRC / f"{k}.npz", allow_pickle=True)
            if abs(float(zz["dr_m"]) - dg) < 0.05:
                keys.append(k)
        sub = sub[sub.key.isin(keys)]
        pick = sub.iloc[rng.choice(len(sub), size=min(N_SAMPLE, len(sub)), replace=False)]
        P, cbh = [], []
        for k, g in pick.groupby("key"):
            zz = np.load(F.SRC / f"{k}.npz", allow_pickle=True)
            kd = np.round(zz["range_rel_m"] / dg).astype(int)
            cols = np.array([int(np.nonzero(kd == kk)[0][0]) for kk in ks])
            idx = np.searchsorted(zz["time_ns"], g.time_ns.to_numpy())
            W = zz["beta_before"][idx][:, cols].astype("float64")
            P.append(W / W.sum(1, keepdims=True))
            cbh.append(np.where(np.isfinite(g.cbh_m), g.cbh_m, 1000.0))
        P, cbh = np.vstack(P), np.concatenate(cbh)
        node = np.clip(np.round((cbh - F.CBH_NODES[0]) / 50.0).astype(int), 0, F.CBH_NODES.size - 1)
        f50o, f100o, kapp_o, _, hh_o = F.observables(P, zg)
        for ns in RESP_NS:
            S = F.peak_aligned(F.gate_sample_conv(conv_gauss(b33, ns), z, dg), ks)   # (alpha, cbh, K)
            S = F._interp_axis(S, np.log(L["alpha_km"]), np.log(F.ALPHA_FINE), axis=0)
            S = F._interp_axis(S, L["cbh_m"], F.CBH_NODES, axis=1)               # (alphaF, cbhN, K)
            S = S / S.sum(-1, keepdims=True)
            ka_m, hh_m = F.model_kapp(np.transpose(S, (1, 0, 2)), zg)           # (cbhN, alphaF)
            al, chi = np.full(P.shape[0], np.nan), np.full(P.shape[0], np.nan)
            for c in np.unique(node):
                m = node == c
                M = S[:, c, :]
                sse = (P[m] ** 2).sum(1)[:, None] + (M ** 2).sum(1)[None, :] - 2 * P[m] @ M.T
                b = sse.argmin(1)
                al[m], chi[m] = F.ALPHA_FINE[b], sse[np.arange(b.size), b]
            # tail-slope alpha: invert the model k_app(alpha) at each profile's CBH node
            ak = np.full(P.shape[0], np.nan)
            for c in np.unique(node):
                m = node == c
                kk = ka_m[c]
                ok = np.isfinite(kk)
                ak[m] = np.interp(kapp_o[m], kk[ok], F.ALPHA_FINE[ok], left=np.nan, right=np.nan)
            rows.append(dict(itype=it, gate_m=dg, response_fwhm_ns=ns, n=P.shape[0],
                             chi2_med=float(np.nanmedian(chi)), alpha_shape_med=float(np.nanmedian(al)),
                             alpha_shape_p10=float(np.nanpercentile(al, 10)),
                             alpha_shape_p90=float(np.nanpercentile(al, 90)),
                             frac_edge=float(np.mean((al <= F.ALPHA_FINE[0]) | (al >= F.ALPHA_FINE[-1]))),
                             alpha_from_kapp_med=float(np.nanmedian(ak)),
                             frac_kapp_outside_grid=float(np.mean(~np.isfinite(ak))),
                             half_height_obs_med=float(np.nanmedian(hh_o)),
                             half_height_model_at_alpha_med=float(np.nanmedian(
                                 hh_m[node, np.abs(F.ALPHA_FINE[None, :] - al[:, None]).argmin(1)]))))
            print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(F.OUT / "response_width_test.csv", index=False, float_format="%.4g")
    print(t.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
