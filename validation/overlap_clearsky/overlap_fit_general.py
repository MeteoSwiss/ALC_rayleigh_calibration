"""General overlap model = sum of K logistic transitions (one may be NEGATIVE = overshoot relaxation),
so a single form reproduces both the monotone factory overlaps and the ones that overshoot >1 then settle.
    O(r) = sum_k  a_k / (1 + exp(-(r - mu_k)/s_k))          O(inf) = sum_k a_k  (~1)
K=2 recovers the double-logistic; K=3-4 adds the overshoot/fine structure."""
import numpy as np
from scipy.optimize import least_squares

def model(r, p):
    K = len(p) // 3; out = np.zeros_like(r, float)
    for k in range(K):
        mu, s, a = p[3*k], p[3*k+1], p[3*k+2]
        out += a / (1.0 + np.exp(-(r - mu) / s))
    return out

def fit_K(r, O, K, zlo, zhi):
    m = np.isfinite(O) & (r >= zlo) & (r <= zhi)
    x = r[m] / 1000.0; y = np.clip(O[m], -0.3, 1.6)          # r in km for conditioning
    # init: spread midpoints across the rise/settle span
    mus = np.linspace(x.min() + 0.05, min(x.max(), 2.2), K)
    p0 = []; lb = []; ub = []
    for k in range(K):
        a0 = 1.0 / K if k < K - 1 else 1.0 - (K - 1) / K
        p0 += [mus[k], 0.12, (0.6 if k == 0 else (a0))]
        lb += [0.0, 0.015, -1.5]; ub += [3.2, 1.5, 2.0]
    def resid(p): return model(x, p) - y
    best = None
    for jit in range(4):                                     # a few restarts with jitter
        pp = np.array(p0, float)
        if jit: pp[0::3] += 0.1 * (jit - 1.5) * np.sign(np.random.default_rng(jit).standard_normal(K)) if False else pp[0::3]
        # deterministic jitter on midpoints
        pp[0::3] = np.clip(np.linspace(x.min()+0.05, min(x.max(),2.0+0.3*jit), K), lb[0], ub[0])
        try:
            r_ = least_squares(resid, pp, bounds=(lb, ub), max_nfev=20000, method='trf')
            rmse = float(np.sqrt(np.mean(r_.fun ** 2)))
            if best is None or rmse < best[1]: best = (r_.x, rmse)
        except Exception:
            continue
    return best  # (params in km-units, rmse)

def eval_km(p, r):  # evaluate a km-fitted param set on a range grid in metres
    return model(r / 1000.0, p)

if __name__ == "__main__":
    import os
    from cfg_tools import parse_cfg
    S = os.path.dirname(__file__)
    cfgs = [("Eriswil", "C:/Users/hervo/OneDrive/Documents/TUB210008_20210630_1024.cfg"),
            ("Melpitz", "C:/Users/hervo/OneDrive/Documents/TUB160061_20161024.cfg"),
            ("W.Canada", "C:/Users/hervo/OneDrive/Documents/overlap/TUB120012_20120917_1024.cfg"),
            ("K.Scheidegg", "C:/Users/hervo/OneDrive/Documents/overlap/TUB120011_20121112_1024.cfg"),
            ("TUB140106", "C:/Users/hervo/OneDrive/Documents/overlap/TUB140106_20150623_1024.cfg")]
    print(f"{'cfg':12s} | {'K=2':>7} {'K=3':>7} {'K=4':>7}   RMSE of sum-of-logistics fit")
    for nm, fp in cfgs:
        v, meta = parse_cfg(fp); res = 15344.0 / len(v); r = np.arange(len(v)) * res
        nz = np.argmax(v > 0.01) * res
        zlo = max(nz - 30, 45); zhi = 2600
        row = []
        for K in [2, 3, 4]:
            fit = fit_K(r, v, K, zlo, zhi); row.append(fit[1] if fit else np.nan)
        print(f"{nm:12s} | {row[0]:7.4f} {row[1]:7.4f} {row[2]:7.4f}")
