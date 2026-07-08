"""Regenerate the figures shown in the top-level README.

    python doc/make_figures.py

Requires the `isodistrreg` Python bindings, numpy, scipy and matplotlib.
All randomness is seeded so the figures are reproducible.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gamma as gamma_dist

from isodistrreg import IDR

# ---------------------------------------------------------------------------
# Palette (a validated, colour-blind-safe categorical set) and shared style.
# ---------------------------------------------------------------------------
BLUE = "#2a78d6"
RED = "#e34948"
AQUA = "#1baf7a"
VIOLET = "#4a3aa7"
ORANGE = "#eb6834"
INK = "#0b0b0b"
MUTED = "#8a8984"

mpl.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        # Computer Modern (the "LaTeX look") via matplotlib's own bundled fonts
        # and mathtext — no external LaTeX/dvipng install required, so the script
        # stays portable. cmr10 has no bold cut, hence normal-weight titles.
        "font.family": "serif",
        "font.serif": ["cmr10", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "axes.formatter.use_mathtext": True,
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.titleweight": "normal",
        "axes.labelcolor": INK,
        "axes.edgecolor": "#d5d4cf",
        "axes.linewidth": 1.0,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": True,
        "grid.color": "#ececea",
        "grid.linewidth": 1.0,
        "legend.frameon": False,
        "figure.dpi": 200,
    }
)


def tidy(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=0)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------
# Figure 1 — a distribution, not a point.
# ---------------------------------------------------------------------------
def figure_distribution():
    rng = np.random.default_rng(20)
    n = 600
    x = rng.uniform(0, 10, n)
    # Mean and spread both grow with x: a genuinely distributional target.
    y = rng.gamma(shape=2.0, scale=(x + 1) / 4)

    fit = IDR(y, x)
    xs = np.linspace(0.2, 9.8, 220)
    qs = fit.quantile(xs[:, None], np.array([0.1, 0.5, 0.9])[None])
    lo, med, hi = qs[:, 0], qs[:, 1], qs[:, 2]
    mean = fit.predict(xs)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 3.9))

    # -- left: scatter + mean + median + central 80% band -----------------
    axL.scatter(x, y, s=9, color=MUTED, alpha=0.35, edgecolors="none",
                label="observations", zorder=1)
    axL.fill_between(xs, lo, hi, color=BLUE, alpha=0.16, linewidth=0,
                     label="10%-90% interval", zorder=2)
    axL.plot(xs, med, color=BLUE, lw=2.2, label="median", zorder=3)
    axL.plot(xs, mean, color=ORANGE, lw=2.2, label="mean", zorder=4)
    axL.set_title(r"One fit $\rightarrow$ the whole conditional distribution")
    axL.set_xlabel(r"covariate $x$")
    axL.set_ylabel(r"response $y$")
    axL.set_ylim(-0.3, min(y.max(), 12))
    axL.legend(loc="upper left")
    tidy(axL)

    # -- right: three predictive CDFs -------------------------------------
    thr = np.asarray(fit.thresholds)
    picks = [2.0, 5.0, 8.0]
    colors = [AQUA, BLUE, VIOLET]
    cdfs = fit.cdf(np.array(picks))
    for xi, c, cdf in zip(picks, colors, cdfs):
        axR.plot(thr, cdf, drawstyle="steps-post", color=c, lw=2.0,
                 label=fr"$x = {xi:.0f}$")
    axR.set_title(r"Predictive CDF $F(y \mid x)$")
    axR.set_xlabel(r"response $y$")
    axR.set_ylabel("probability")
    axR.set_xlim(0, 10)
    axR.set_ylim(-0.02, 1.02)
    axR.legend(loc="lower right", title="fitted at")
    tidy(axR)

    fig.tight_layout()
    fig.savefig("doc/idr_distribution.png", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 — right-censoring.
# ---------------------------------------------------------------------------
def figure_censoring():
    rng = np.random.default_rng(4)
    n = 1500
    x = rng.uniform(0, 10, n)
    event = rng.gamma(shape=2.0, scale=(x + 1) / 4)       # true event time
    cens = rng.exponential(scale=5.0, size=n)             # censoring time
    t = np.minimum(event, cens)
    observed = event <= cens                              # True = event seen

    # Survival IDR uses the censoring indicator; the naive fit ignores it and
    # treats every censored time as if the event happened right then.
    fit_sidr = IDR(t, x, observed)
    fit_naive = IDR(t, x)

    x0 = 6.0
    thr_s = np.asarray(fit_sidr.thresholds)
    thr_n = np.asarray(fit_naive.thresholds)
    surv_sidr = 1 - fit_sidr.cdf(np.array([x0]))[0]
    surv_naive = 1 - fit_naive.cdf(np.array([x0]))[0]
    grid = np.linspace(0, thr_n.max(), 400)
    surv_true = gamma_dist.sf(grid, a=2.0, scale=(x0 + 1) / 4)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(grid, surv_true, color=INK, lw=1.6, ls=(0, (4, 3)),
            label="true survival")
    ax.plot(thr_n, surv_naive, drawstyle="steps-post", color=RED, lw=2.2,
            label="naive (censoring ignored)")
    ax.plot(thr_s, surv_sidr, drawstyle="steps-post", color=BLUE, lw=2.2,
            label="S-IDR (censoring handled)")
    ax.set_title(fr"Survival estimate at $x = {x0:.0f}$   ({100*(1-observed.mean()):.0f}% censored)")
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"P(event after $t \mid x$)")
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="upper right")
    tidy(ax)

    fig.tight_layout()
    fig.savefig("doc/idr_censoring.png", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3 — subagging.
# ---------------------------------------------------------------------------
def figure_subagging():
    rng = np.random.default_rng(11)
    n = 120
    x = rng.uniform(0, 10, n)
    y = rng.gamma(shape=2.0, scale=(x + 1) / 4)
    raw = IDR(y, x)
    bag = IDR(y, x, subsamples=50, subsample_size=0.5, seed=1, n_jobs=4)

    xs = np.linspace(0.2, 9.8, 300)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(xs, raw.predict(xs), color=MUTED, lw=1.8, label="single fit")
    ax.plot(xs, bag.predict(xs), color=BLUE, lw=2.4, label=r"subagged ($50\!\times\!$)")
    ax.set_title("Subagging stabilises the fit")
    ax.set_xlabel(r"covariate $x$")
    ax.set_ylabel("predicted mean")
    ax.set_xlim(0, 10)
    ax.legend(loc="upper left")
    tidy(ax)

    fig.tight_layout()
    fig.savefig("doc/idr_subagging.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    figure_distribution()
    figure_censoring()
    figure_subagging()
    print("wrote doc/idr_distribution.png, idr_censoring.png, idr_subagging.png")
