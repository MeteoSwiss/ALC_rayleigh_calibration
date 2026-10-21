# -*- coding: utf-8 -*-
"""Illustrations for the CHM15k temperature-calibration deck (report 19). Real numbers only.

Run with the project venv:  .venv/Scripts/python.exe presentation/chm15k_temperature/make_illustrations.py
"""
import sys
from pathlib import Path

import numpy as np

SKILL = Path("C:/Users/hervo/.claude/plugins/cache/presentation/meteoswiss-presentation/0.1.0/"
             "skills/meteoswiss-presentation/scripts")
sys.path.insert(0, str(SKILL))
from illustration import C, blank, box, arrow, new_figure, save  # noqa: E402

HERE = Path(__file__).resolve().parent
FIG = HERE / "fig"


def timeline():
    """Slide 2: night-only calibration vs the daytime warming and the resulting excess."""
    z = np.load(FIG / "diurnal_composite.npz")
    h = np.arange(24) + 0.5
    fig, ax = new_figure("figure")
    ax.axvspan(0, 4, color=C["blue_light"], lw=0)
    ax.axvspan(22, 24, color=C["blue_light"], lw=0)
    ax.text(2.0, 30.2, "Rayleigh\ncalibration\n(night only)", ha="center", va="top", fontsize=12, color=C["blue_dark"])
    ax.text(23.0, 30.2, "night", ha="center", va="top", fontsize=12, color=C["blue_dark"])
    ax.plot(h, z["T_jja"], "-", color=C["brown"], lw=2.5, label="internal temperature rise [K], summer median")
    ax.plot(h, z["E_jja"], "-", color=C["red"], lw=3, label="excess of the effective constant [%], summer median")
    ax.plot(h, z["E_p90_jja"], "--", color=C["red"], lw=1.6, label="same, 90th percentile of units")
    ax.set_xlim(0, 24)
    ax.set_ylim(-1, 31)
    ax.set_xticks(range(0, 25, 3))
    ax.set_xlabel("local solar hour")
    ax.set_ylabel("relative to night")
    ax.legend(loc="upper left", bbox_to_anchor=(0.19, 1.0), frameon=False,
              framealpha=0.9, fontsize=11)
    save(fig, FIG / "fig_timeline.png")


def estimators():
    """Slide 4: the temperature coefficient from independent estimators (real values, 1-sigma)."""
    rows = [  # label, value, err (or (lo, hi) range), colour
        ("Naive regression (season mixed in)", 2.78, 0.03, C["grey"]),
        ("Warm/cold spells, 158 units 2025-26", 2.06, 0.04, C["blue_dark"]),
        ("  + aerosol control (AERONET AOD, 27 units)", 2.01, 0.11, C["blue_dark"]),
        ("Warm/cold spells, 11-yr archive, 13 stations", 1.43, 0.09, C["blue_dark"]),
        ("Within single nights, 150 units", 1.46, 0.02, C["blue_dark"]),
        ("Co-located pairs at night (atmosphere cancels)", (0.4, 2.1), None, C["green_light"]),
    ]
    fig, ax = new_figure("figure")
    y = np.arange(len(rows))[::-1]
    for yi, (lab, v, e, col) in zip(y, rows):
        if isinstance(v, tuple):
            ax.plot(v, [yi, yi], "-", color=col, lw=10, solid_capstyle="butt")
            ax.text(v[1] + 0.08, yi, f"+{v[0]:.1f} ... +{v[1]:.1f}", va="center", fontsize=12)
        else:
            ax.barh(yi, v, height=0.55, color=col)
            ax.errorbar(v, yi, xerr=e, fmt="none", ecolor="black", capsize=3, lw=1)
            ax.text(v + e + 0.06, yi, f"+{v:.2f}", va="center", fontsize=12)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=12)
    ax.set_xlim(0, 3.4)
    ax.set_xlabel("temperature coefficient κ of the CHM15k constant [%/K]")
    ax.axvline(0, color=C["grey"], lw=0.8)
    save(fig, FIG / "fig_estimators.png")


def mechanism():
    """Slide 9: candidate causes and verdicts (text_image slot)."""
    rows = [
        ("Aerosol in the calibration window", "excluded", "AERONET: 2.36 → 2.01 %/K"),
        ("Detector / reference pulse", "excluded", "raw counts keep the slope"),
        ("Electronic baseline (dark)", "excluded", "same on clean & hazy nights"),
        ("Overlap / optical alignment", "excluded", "flat from 0.5 to 8 km"),
        ("Laser pulse energy", "unlikely", "pulse rate −0.06 %/K"),
        ("Laser line vs narrow filter", "favoured", "curved response, 0.1 nm filter"),
    ]
    col = {"excluded": C["grey_light"], "unlikely": C["yellow_light"], "favoured": C["red"]}
    txt = {"excluded": "black", "unlikely": "black", "favoured": "white"}
    fig, ax = new_figure("text_image")
    blank(ax)
    n = len(rows)
    for k, (cause, verdict, why) in enumerate(rows):
        yc = 1 - (k + 0.5) / n
        ax.text(0.0, yc + 0.035, cause, ha="left", va="center", fontsize=12, fontweight="bold")
        ax.text(0.0, yc - 0.045, why, ha="left", va="center", fontsize=10.5, color=C["grey"])
        box(ax, 0.86, yc, verdict, w=0.26, h=0.11, fc=col[verdict], fontsize=11.5, color=txt[verdict],
            bold=verdict == "favoured")
    save(fig, FIG / "fig_mechanism.png")


def roadmap():
    """Slide 10: staged introduction of the correction (text_image slot)."""
    fig, ax = new_figure("text_image")
    blank(ax)
    steps = [
        (0.83, "1  Calibration side (now)", "κ per module, publish C_ref and κ\nno product change", C["blue_light"]),
        (0.53, "2  Shadow L2 (3–6 months)", "pair, CL61 and EARLINET sites\ngo / no-go on day−night bias", C["yellow_light"]),
        (0.23, "3  L1→L2 hub, unit by unit", "only units with a reliable κ\nκ = 0 stays the default", C["green_light"]),
    ]
    for y, head, body, fc in steps:
        box(ax, 0.5, y, "", w=0.98, h=0.24, fc=fc)
        ax.text(0.04, y + 0.06, head, ha="left", va="center", fontsize=13, fontweight="bold")
        ax.text(0.04, y - 0.05, body, ha="left", va="center", fontsize=11)
    for y0, y1 in ((0.71, 0.65), (0.41, 0.35)):
        arrow(ax, (0.5, y0), (0.5, y1))
    ax.text(0.5, 0.03, "in parallel: questions to Lufft (filter, laser energy, module generations)",
            ha="center", va="center", fontsize=10.5, color=C["brown"])
    save(fig, FIG / "fig_roadmap.png")


def filter_line():
    """Filter passband drifting against a fixed laser line (schematic, text_image slot)."""
    import matplotlib.pyplot as plt
    fig0, _ = new_figure("text_image")          # sets the house rcParams
    plt.close(fig0)
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(4.55, 3.3), gridspec_kw=dict(height_ratios=[1.7, 1]))
    x = np.linspace(-2.2, 2.2, 400)

    def passband(c):
        # flat-topped (super-Gaussian) interference-filter shape, schematic
        return np.exp(-((x - c) / 0.8) ** 6)
    line = 0.45
    ax.plot(x, passband(-0.35), "-", color=C["blue_dark"], lw=2.2, label="filter, cool")
    ax.plot(x, passband(0.0), "-", color=C["red"], lw=2.2, label="filter, warm")
    ax.axvline(line, color="black", lw=2.2)
    ax.text(line + 0.08, 1.33, "laser line\n(regulated, fixed)", fontsize=10, va="top")
    ax.annotate("", xy=(-0.1, 0.5), xytext=(-0.9, 0.5),
                arrowprops=dict(arrowstyle="-|>", color=C["grey"], lw=1.6))
    ax.text(-1.72, 0.5, "warmer:\npassband\nshifts →", fontsize=9.5, ha="center", va="center",
            color=C["brown"])
    ax.plot([line], [np.interp(line, x, passband(-0.35))], "o", color=C["blue_dark"], ms=6)
    ax.plot([line], [np.interp(line, x, passband(0.0))], "o", color=C["red"], ms=6)
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(0, 1.35)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("wavelength (schematic; datasheet passband 1 nm)", fontsize=10)
    ax.set_ylabel("transmission", fontsize=10)
    ax.legend(frameon=False, fontsize=9.5, loc="upper left", ncol=2, bbox_to_anchor=(-0.02, 1.04),
              handlelength=1.2, columnspacing=0.8)
    tt = np.linspace(-0.6, 0.6, 50)
    tr = np.array([np.exp(-((line - c) / 0.8) ** 6) for c in tt])
    ax2.plot(tt, tr, "-", color=C["red"], lw=2.2)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xlabel("instrument temperature", fontsize=10)
    ax2.set_ylabel("signal", fontsize=10)
    ax2.text(0.98, 0.08, "rises, then flattens\n(curvature seen in the data)", transform=ax2.transAxes,
             fontsize=9.5, ha="right", va="bottom", color=C["brown"])
    save(fig, FIG / "fig_filter_line.png")


def ml_vs_exp():
    """New slide: prediction error of the simple exponential vs the housekeeping ML (2024 test year)."""
    import csv
    r = {row["model"]: float(row["robust_sd_common"]) for row in
         csv.DictReader(open("C:/DATA/Projects/202610_CHM15k_stability/ml_vs_exp_scores.csv", encoding="utf-8"))}
    groups = [("housekeeping only", [("ML\n(bagged tree)", r["ML_A"], C["grey"]),
                                     ("exp(κΔT)", r["EXP_A"], C["red"])]),
              ("with last 30 days\nof calibrations", [("running\nmedian", r["RM_B"], C["grey_light"]),
                                                     ("ML ×\nbias", r["ML_B"], C["grey"]),
                                                     ("+ exp(κΔT)", r["EXP_B"], C["red"])])]
    fig, ax = new_figure("text_image")
    xpos, labels, x = [], [], 0
    for g, items in groups:
        for lab, v, col in items:
            ax.bar(x, v, width=0.75, color=col, edgecolor=C["grey"], lw=0.6)
            ax.text(x, v + 0.4, f"{v:.1f}", ha="center", va="bottom", fontsize=11)
            xpos.append(x)
            labels.append(lab)
            x += 1
        ax.text(x - len(items) / 2 - 0.5, 27.5, g, ha="center", va="top", fontsize=10.5, color=C["brown"])
        x += 0.7
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylim(0, 28)
    ax.set_ylabel("prediction error on 2024 [%]", fontsize=11)
    save(fig, FIG / "fig_ml_vs_exp.png")


def _math_rc():
    """STIX sans for mathtext: close to Arial, and it has scalable delimiters (Arial has none)."""
    import matplotlib.pyplot as plt
    plt.rcParams.update({"mathtext.fontset": "stixsans"})


def equation_vs_ml():
    """Slide 10: the correction equation (typeset) next to the ML comparison (figure slot)."""
    import csv
    import matplotlib.pyplot as plt
    fig0, _ = new_figure("figure")
    plt.close(fig0)
    _math_rc()
    fig = plt.figure(figsize=(8.8, 2.95))
    axe = fig.add_axes([0.0, 0.0, 0.47, 1.0])
    blank(axe)
    axe.text(0.02, 0.86, r"$C_L(t) = C_{\mathrm{ref}}\;\exp\!\left[\kappa_u\,\left(T_{\mathrm{int}}(t) - T_{\mathrm{ref}}\right)\right]$",
             fontsize=17, va="center")
    lines = [(r"$\kappa_u$", "temperature coefficient of the unit [%/K],\none per optical module, from clear nights"),
             (r"$T_{\mathrm{int}}(t)$", "internal temperature, already in L1"),
             (r"$C_{\mathrm{ref}}$", "nightly Rayleigh constant normalised\nto the reference temperature"),
             (r"$T_{\mathrm{ref}}$", "fixed reference temperature (e.g. 25 °C)")]
    y = 0.64
    for sym, txt in lines:
        axe.text(0.04, y, sym, fontsize=13, va="top")
        axe.text(0.20, y, txt, fontsize=11, va="top", color=C["black"])
        y -= 0.165 if "\n" in txt else 0.11
    r = {row["model"]: float(row["robust_sd_common"]) for row in
         csv.DictReader(open("C:/DATA/Projects/202610_CHM15k_stability/ml_vs_exp_scores.csv", encoding="utf-8"))}
    ax = fig.add_axes([0.56, 0.22, 0.43, 0.68])
    items = [("ML", r["ML_A"], C["grey"], 0), ("exp", r["EXP_A"], C["red"], 1),
             ("running\nmedian", r["RM_B"], C["grey_light"], 2.7), ("ML ×\nbias", r["ML_B"], C["grey"], 3.7),
             ("+ exp", r["EXP_B"], C["red"], 4.7)]
    for lab, v, col, x in items:
        ax.bar(x, v, width=0.8, color=col, edgecolor=C["grey"], lw=0.6)
        ax.text(x, v + 0.5, f"{v:.1f}", ha="center", va="bottom", fontsize=10)
    ax.set_xticks([it[3] for it in items])
    ax.set_xticklabels([it[0] for it in items], fontsize=9.5)
    ax.text(0.5, 27.5, "housekeeping only", ha="center", fontsize=9.5, color=C["brown"])
    ax.text(3.7, 27.5, "with last 30 days of calibrations", ha="center", fontsize=9.5, color=C["brown"])
    ax.set_ylim(0, 29)
    ax.set_ylabel("prediction error 2024 [%]", fontsize=10.5)
    ax.tick_params(axis="y", labelsize=10)
    fig.savefig(FIG / "fig_equation_ml.png", dpi=200, facecolor="white")
    plt.close(fig)


def model_slide():
    """Model slide: equations, inputs with sources and uncertainties, and the predicted response band."""
    import csv
    import matplotlib.pyplot as plt
    fig0, _ = new_figure("figure")
    plt.close(fig0)
    _math_rc()
    rows = [r for r in csv.DictReader(open("C:/DATA/Projects/202610_CHM15k_stability/q11_filter_inversion.csv",
                                             encoding="utf-8")) if r.get("alpha_over_W") not in (None, "")]
    aw = np.array([float(r["alpha_over_W"]) for r in rows])
    d0 = np.array([float(r["delta0_over_W"]) for r in rows])
    tm = float(np.median([float(r["T_mean"]) for r in rows]))
    fig = plt.figure(figsize=(8.8, 3.3))
    axe = fig.add_axes([0.0, 0.80, 1.0, 0.20])
    blank(axe)
    axe.text(0.0, 0.55, r"$S\ \propto\ \mathcal{T}(\lambda_L)=\exp\left[-4\ln 2\,\left(\frac{\lambda_L-\lambda_c(T_f)}{W}\right)^{2}\right]$",
             fontsize=13, va="center")
    axe.text(0.56, 0.55, r"$\lambda_c(T_f)=\lambda_{c0}+\alpha\,(T_f-T_0)\qquad T_f\approx T_{\mathrm{int}}$",
             fontsize=13, va="center")
    ax = fig.add_axes([0.065, 0.14, 0.31, 0.55])
    T = np.linspace(10, 40, 121)
    curves = np.array([100 * (np.exp(-4 * np.log(2) * ((d - a * (T - tm)) ** 2 - d ** 2)) - 1)
                       for a, d in zip(aw, d0)])
    ax.axvspan(15, 26, color=C["grey_light"], lw=0)
    ax.fill_between(T, np.percentile(curves, 25, axis=0), np.percentile(curves, 75, axis=0), color=C["red"],
                    alpha=0.2, lw=0, label="IQR of 110 units")
    ax.plot(T, np.median(curves, axis=0), color=C["red"], lw=2, label="median unit")
    ax.axvline(33.85, color=C["brown"], ls="--", lw=1.2)
    ax.text(33.4, -27, "module\nset-point", fontsize=8, color=C["brown"], ha="right", va="bottom")
    ax.text(20.5, -27, "operating range", fontsize=8, ha="center", va="bottom", color=C["grey"])
    ax.axhline(0, color=C["grey"], lw=0.6)
    ax.set_xlim(10, 40)
    ax.set_ylim(-28, 22)
    ax.set_xlabel("filter temperature [°C]", fontsize=9.5)
    ax.set_ylabel("signal vs mean [%]", fontsize=9.5)
    ax.tick_params(labelsize=8.5)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    axt = fig.add_axes([0.42, 0.0, 0.58, 0.78])
    blank(axt)
    head = ["input", "value", "source / uncertainty"]
    data = [
        (r"filter FWHM $W$", "1 nm (assumed)", "Lufft datasheet 2024 (manual: '1 nm or less')"),
        (r"filter drift $\alpha$", "2–30 pm/K", "hard 2–5; conventional 14–30 (suppliers)"),
        (r"laser line $\lambda_L$", "fixed", "LOM heated (manual R19); +3–5 pm/K if leaking"),
        (r"coupling $T_f$ / $T_{\mathrm{int}}$", "1, lag < 0.5 h", "this study; filter is inside the LOM?"),
        ("passband shape", "Gaussian", "assumed; a flat-top changes the curvature"),
        (r"fit per unit: $\kappa$, curvature", "150 / 141 units", "curvature significant in only 9"),
        (r"result: $\alpha$ at $W$ = 1 nm", "≈ 20 pm/K", "IQR 15–24; Jenoptik filter spec: 20 pm/K"),
        (r"result: detuning $\Delta_0/W$", "+0.19", "IQR 0.11–0.32; ~10 % transmission loss"),
        (r"result: best-match $T^{*}$", "≈ 31 °C", "IQR 26–43 °C; module set-point 33.85 °C"),
    ]
    xs = [0.0, 0.29, 0.495]
    y = 0.96
    for x, h in zip(xs, head):
        axt.text(x, y, h, fontsize=9, fontweight="bold", va="center")
    axt.plot([0, 1], [y - 0.05, y - 0.05], color=C["grey"], lw=0.8)
    for row in data:
        y -= 0.103
        for x, cell in zip(xs, row):
            axt.text(x, y, cell, fontsize=8.3, va="center")
    fig.savefig(FIG / "fig_model.png", dpi=200, facecolor="white")
    plt.close(fig)


def filters_vs_fit():
    """Commercial 1064 nm filters in the (width, thermal drift) plane against the band allowed by the fit."""
    from matplotlib.patches import Rectangle
    _math_rc()
    fig, ax = new_figure("figure")
    W = np.logspace(-1, np.log10(20), 50)
    # observed alpha/W = 0.015-0.024 /K (q11 inversion IQR); alpha in pm/K, W in nm
    ax.fill_between(W, 15 * W, 24 * W, color=C["red"], alpha=0.18, lw=0)
    ax.plot(W, 20 * W, color=C["red"], lw=2)
    ax.text(0.108, 14, "fit to the CHM15k data\n(κ ≈ +1.5…+2 %/K)", color=C["red"], fontsize=11.5, va="center")
    for k in (0.2, 0.5):
        ax.plot(W, k * 9.49 * W, ":", color=C["grey"], lw=1)
        ax.text(3.2, k * 9.49 * 3.2 * 1.15, f"κ = {k:g} %/K", color=C["grey"], fontsize=9.5, ha="center",
                va="bottom", rotation=13)
    ax.axvline(1.0, color=C["grey"], ls="--", lw=1)
    ax.text(1.06, 1.12, "Lufft datasheet: 1 nm (manual: '1 nm or less')", fontsize=9.5, color=C["grey"])
    ax.errorbar([1.0], [30], xerr=[[0.2], [0.2]], yerr=[[2], [2]], fmt="s", color=C["brown"], ms=9, capsize=3)
    ax.text(0.76, 42, "Andover standard 1 nm\n(~30 pm/K, generic chart)", fontsize=10.5, color=C["brown"],
            ha="right", va="center")
    ax.add_patch(Rectangle((0.4, 2), 0.6, 3, facecolor=C["blue_light"], edgecolor=C["blue_dark"], lw=1.5))
    ax.text(0.63, 3.15, "Alluxa ultra-narrow\nhard-coated", fontsize=10, color=C["blue_dark"], ha="center",
            va="center", zorder=5, bbox=dict(facecolor=C["blue_light"], edgecolor="none", pad=1))
    ax.plot([10, 10], [2, 5], color=C["blue_dark"], lw=5, solid_capstyle="butt")
    ax.text(9.2, 3.2, "hard-coated\n10 nm", fontsize=10.5, color=C["blue_dark"], ha="right", va="center")
    ax.annotate("", xy=(10, 9), xytext=(10, 20), arrowprops=dict(arrowstyle="->", color=C["grey"], lw=1.5))
    ax.plot([10], [20], "_", color=C["grey"], ms=12, mew=2)
    ax.text(11, 12, "Newport 10 nm\n< 20 pm/K", fontsize=9.5, color=C["grey"], va="center")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 20)
    ax.set_ylim(1, 60)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10, 20])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10", "20"])
    ax.set_yticks([1, 2, 5, 10, 20, 50])
    ax.set_yticklabels(["1", "2", "5", "10", "20", "50"])
    ax.set_xlabel("filter width FWHM [nm]")
    ax.set_ylabel("filter drift dλ/dT [pm/K]")
    save(fig, FIG / "fig_filters_vs_fit.png")


def filters_alpha():
    """With W = 1 nm: filter drift required by the fit vs the drift of commercial filter technologies."""
    _math_rc()
    fig, ax = new_figure("figure")
    rows = [  # label, low, high, kind ("pt", "range", "max"), colour
        ("Jenoptik JENOFILT (Jena), spec", 20, 20, "pt", C["brown"]),
        ("Andover standard bandpass, chart", 28, 32, "range", C["brown"]),
        ("Blifford 1966, commercial, measured", 8, 29, "range", C["brown"]),
        ("Omega, unassisted coating (VIS)", 14, 14, "pt", C["brown"]),
        ("Newport 1064 / Thorlabs FL, spec", 0, 20, "max", C["grey"]),
        ("Omega, plasma-sputtered (VIS)", 7, 7, "pt", C["blue_dark"]),
        ("Semrock MaxLine 1064, spec", 0, 5.3, "max", C["blue_dark"]),
        ("Alluxa hard-coated, spec", 2, 5, "range", C["blue_dark"]),
        ("Tsai 2001, sputtered, measured", 0, 3, "max", C["blue_dark"]),
    ]
    ax.axvspan(15, 24, color=C["red"], alpha=0.18, lw=0)
    ax.axvline(20, color=C["red"], lw=2)
    ax.text(20.4, len(rows) - 0.35, "required by the fit\nfor W = 1 nm", color=C["red"], fontsize=11, va="top")
    for i, (lab, lo, hi, kind, col) in enumerate(rows):
        y = len(rows) - 1 - i
        if kind == "pt":
            ax.plot([lo], [y], "o", color=col, ms=9)
        elif kind == "range":
            ax.plot([lo, hi], [y, y], color=col, lw=6, solid_capstyle="butt")
        else:
            ax.annotate("", xy=(0.3, y), xytext=(hi, y), arrowprops=dict(arrowstyle="->", color=col, lw=1.8))
            ax.plot([hi], [y], "|", color=col, ms=12, mew=2.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=10.5)
    for t, r in zip(ax.get_yticklabels(), rows[::-1]):
        t.set_color(r[4])
    ax.set_xlim(0, 35)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("filter drift dλ/dT [pm/K]")
    top = ax.secondary_xaxis("top", functions=(lambda a: 0.1054 * a, lambda k: k / 0.1054))
    top.set_xlabel("predicted κ for a 1 nm filter [%/K]", fontsize=11)
    top.tick_params(labelsize=10)
    ax.text(34.5, 6.0, "conventional\n(cemented / evaporated)", color=C["brown"], fontsize=10, ha="right")
    ax.text(34.5, 1.0, "hard-coated\n(sputtered)", color=C["blue_dark"], fontsize=10, ha="right")
    save(fig, FIG / "fig_filters_alpha.png")




if __name__ == "__main__":
    filters_alpha()
    filters_vs_fit()
    equation_vs_ml()
    model_slide()
    timeline()
    estimators()
    mechanism()
    roadmap()
    filter_line()
    ml_vs_exp()
    print("ok")
