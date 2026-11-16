"""Slide 9 (Annual Meeting 2026): C_L vs cloud-base height at Payerne, 4 eta/WV configurations.
Reads presentation/figs_annual2026/eta_ladder/*.npz (cloud_profile_dump outputs).
Usage: python presentation/slide9_eta_ladder_4cfg.py [out.png]"""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = sys.argv[2] if len(sys.argv) > 2 else 'presentation/figs_annual2026/eta_ladder/'
OUT = sys.argv[1] if len(sys.argv) > 1 else 'presentation/figs_annual2026/fig_eta_ladder_4cfg.png'

def load(f):
    z = np.load(D + f + '.npz', allow_pickle=True)
    cbh = z['cbh_m'] / 1000.; ci = 1 / z['c_oconnor']; d = np.asarray(z['date'])
    ok = np.isfinite(cbh) & np.isfinite(ci) & (ci > 0)
    return cbh[ok], ci[ok], d[ok]

def cl_ols(x, y, g):  # OLS slope with day-clustered (CR1) standard error
    X = np.c_[np.ones_like(x), x]; b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b
    XtXi = np.linalg.inv(X.T @ X); G = np.unique(g); meat = np.zeros((2, 2))
    for gg in G:
        m = g == gg; s = X[m].T @ r[m]; meat += np.outer(s, s)
    n, k, ng = len(y), 2, len(G); c = ng / (ng - 1) * (n - 1) / (n - k)
    V = c * XtXi @ meat @ XtXi
    return b[1], np.sqrt(V[1, 1])

cfg = [('hopkin_nowv', 'Hopkin η, no WV', '#A6A6A6', 's'),
       ('hopkin_wv',   'Hopkin η + WV',   '#303030', 'v'),
       ('nowv',        'New η, no WV',    '#56B4E9', 'o'),
       (None,          'New η + WV (oper.)', '#D55E00', 'D')]
inst = [('B', 'CL31 Payerne', 'nominal', '909.7 nm / 6 nm'),
        ('C', 'CL61 Payerne', 'l55s008', '910.55 nm / 0.19 nm')]
edges = np.arange(0.5, 2.51, 0.25)
fig, axs = plt.subplots(1, 2, figsize=(17.28, 5.9), sharey=True, constrained_layout=True)
for ax, (ident, name, optag, line) in zip(axs, inst):
    for tag, lab, col, mk in cfg:
        t = tag or optag
        cbh, ci, d = load(f'{t}_06610_{ident}')
        pct = 100 * ci / np.median(ci); w = (pct >= 60) & (pct <= 140)
        s, se = cl_ols(cbh[w], 100 * np.log(ci[w]), d[w])
        yc, m, e = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = w & (cbh >= lo) & (cbh < hi); v = pct[sel]
            if v.size >= 20:
                yc.append((lo + hi) / 2); m.append(v.mean())
                e.append(v.std() / np.sqrt(max(1, len(np.unique(d[sel])))))
        yc, m, e = map(np.array, (yc, m, e))
        ax.fill_betweenx(yc, m - e, m + e, color=col, alpha=0.13, lw=0)
        ax.plot(m, yc, '-', marker=mk, color=col, lw=2.6 if tag is None else 2, ms=8,
                label=f'{lab}   {s:+.1f} ± {se:.1f}'.replace('-', '−'))
        print(name, t, round(s, 2), round(se, 2), len(np.unique(d[w])), int(w.sum()))
    ax.axvline(100, color='k', ls='--', lw=0.9)
    ax.set_xlim(88, 110); ax.set_ylim(0.4, 2.5)
    ax.set_title(f'{name}  (WV line {line})', fontsize=19)
    ax.set_xlabel('$C_L$ / unit median  [%]', fontsize=17)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), fontsize=15, frameon=False, ncol=2,
              title='slope of $C_L$ with cloud base  [%/km]', title_fontsize=15, columnspacing=1.6)
    ax.grid(alpha=0.25); ax.tick_params(labelsize=15)
axs[0].set_ylabel('Cloud-base height  [km]', fontsize=17)
fig.savefig(OUT, dpi=150)
print('saved', OUT)
