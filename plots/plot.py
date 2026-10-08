"""
Written by Claude
The figures for the write-up, as PDF and PNG in figures/.

1. z_curve: Z(t) with its first zeros marked (data/z_curve.csv, from the production z(t); zeros
   from the zero finder's own output).
2. term_errors: error of Z(t) against the number of correction terms, against mpmath
   (data/term_errors_summary.csv).
3. count_staircase: the number of zeros up to T (a step function from the zero list) against the
   smooth Riemann-von Mangoldt count, over the whole zero list, with an inset for T up to 100.
4. spacings: gaps between consecutive zeros, normalised by the local mean spacing, against the
   Wigner surmise for the GUE, p(s) = (32 / pi^2) s^2 exp(-4 s^2 / pi) (Wikipedia, "Wigner surmise";
   checked: it integrates to 1 with mean 1).

Reads data/zeros.csv (the main program's output, e.g. from `make run N=1000000`) and the files
`make plots` writes, which also calls this script.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

# Colours: the first three slots of the reference palette, validated for all-pairs colour-blind
# separation on a light surface (aqua is below 3:1 contrast, hence direct labels on every line).
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_SECONDARY, GRID = "#0b0b0b", "#52514e", "#e6e5e1"

ZEROS_CSV = "data/zeros.csv"
OUT_DIR = "figures"

plt.rcParams.update({
    "figure.figsize": (6.4, 3.6),
    "font.size": 10,
    "text.color": INK,
    "axes.labelcolor": INK_SECONDARY,
    "axes.edgecolor": INK_SECONDARY,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.axisbelow": True,
    "legend.frameon": False,
    "savefig.bbox": "tight",
})


def e_notation(value, _position=None):
    """Tick label like 2e5 or 1e6, written on every tick (plain text, so it stays narrow)."""
    if value == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = round(value / 10**exponent, 6)
    return f"{mantissa:g}e{exponent}"


def save(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(f"{OUT_DIR}/{name}.pdf")
    fig.savefig(f"{OUT_DIR}/{name}.png", dpi=600)
    plt.close(fig)
    print(f"wrote {OUT_DIR}/{name}.pdf and .png")


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def figure_z_curve():
    curve = read_csv("data/z_curve.csv")
    zeros = read_csv(ZEROS_CSV)["t"]
    zeros = zeros[zeros <= curve["t"].max()]

    fig, ax = plt.subplots()
    ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    ax.plot(curve["t"], curve["z"], color=BLUE, linewidth=1.6)
    # Zeros sit on the y = 0 line; a 2 px surface ring keeps the markers readable over the curve.
    ax.plot(zeros, np.zeros_like(zeros), "o", color=INK, markersize=5,
            markeredgecolor="white", markeredgewidth=1.0, zorder=3)
    ax.annotate(f"first zero, t = {zeros[0]:.3f}", xy=(zeros[0], 0), xytext=(zeros[0] + 1.5, 0.7 * curve["z"].min()),
                color=INK_SECONDARY, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK_SECONDARY, linewidth=0.8))
    ax.set_xlabel("t")
    ax.set_ylabel("Z(t)")
    ax.set_title(f"Z(t) and its first {len(zeros)} zeros", loc="left", fontsize=10, color=INK)
    ax.set_xlim(curve["t"].min(), curve["t"].max())
    save(fig, "z_curve")


def figure_term_errors():
    data = read_csv("data/term_errors_summary.csv")
    t = data["t"]
    series = [
        ("c0", r"$C_0$", BLUE, -0.75, r"$t^{-3/4}$"),
        ("c01", r"$C_0$–$C_1$", ORANGE, -1.25, r"$t^{-5/4}$"),
        ("c012", r"$C_0$–$C_2$", AQUA, -1.75, r"$t^{-7/4}$"),
    ]

    fig, ax = plt.subplots()
    for key, label, colour, slope, slope_label in series:
        err = data[key]
        # Reference line: the predicted power law, from the first point of this series.
        ax.plot(t, err[0] * (t / t[0]) ** slope, "--", color=INK_SECONDARY, linewidth=0.8, zorder=1)
        ax.plot(t, err, "o-", color=colour, linewidth=1.6, markersize=5,
                markeredgecolor="white", markeredgewidth=1.0, label=label, zorder=2)
        ax.annotate(label, xy=(t[-1], err[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", color=INK, fontsize=9)
        ax.annotate(slope_label, xy=(t[2], err[0] * (t[2] / t[0]) ** slope),
                    xytext=(-4, -12), textcoords="offset points", ha="right",
                    color=INK_SECONDARY, fontsize=8)

    # The C0-C2 error turns back up once double-precision rounding of the phase exceeds it. The note
    # sits in the wedge between the rising C0-C2 line and the falling dashed power law, which widens
    # to the right. No leader line: the text is next to the point it describes.
    ax.text(2.8e5, 1.25e-10, "double-precision\nrounding takes over", ha="center", va="center",
            color=INK_SECONDARY, fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("t")
    ax.set_ylabel("Max |error| in Z(t), over p")
    ax.set_title("Error against number of correction terms", loc="left", fontsize=10, color=INK)
    ax.set_xlim(t[0] * 0.9, t[-1] * 4)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], linestyle="--", color=INK_SECONDARY, linewidth=0.8))
    labels.append("Predicted power law")
    ax.legend(handles, labels, loc="lower left")
    save(fig, "term_errors")


def smooth_count(T):
    """Riemann-von Mangoldt: (T / 2 pi) ln(T / 2 pi) - T / 2 pi + 7/8."""
    x = T / (2.0 * np.pi)
    return x * np.log(x) - x + 7.0 / 8.0


def figure_count_staircase():
    zeros = read_csv(ZEROS_CSV)["t"]
    n = np.arange(1, len(zeros) + 1)

    fig, ax = plt.subplots()
    # The two curves coincide at this scale, so the formula is drawn thick underneath and the zero
    # count thin on top. Plotting every 1000th zero keeps the file small; the steps are invisible anyway.
    grid = np.linspace(zeros[0], zeros[-1], 2000)
    ax.plot(grid, smooth_count(grid), color=ORANGE, linewidth=4.0, solid_capstyle="butt",
            label="RvM formula")
    ax.plot(zeros[::1000], n[::1000], color=BLUE, linewidth=1.2, label="Zeros found")
    ax.legend(loc="lower right")
    ax.set_xlabel("T")
    ax.set_ylabel("number of zeros up to T")
    ax.set_title(f"All zeros found, against the Riemann-von Mangoldt formula",
                 loc="left", fontsize=10, color=INK)
    ax.set_xlim(0, zeros[-1])
    ax.set_ylim(0, None)
    ax.xaxis.set_major_formatter(FuncFormatter(e_notation))
    ax.yaxis.set_major_formatter(FuncFormatter(e_notation))

    # Inset: the same two curves for T up to 100, where the steps are visible.
    t_max = 100.0
    inside = zeros <= t_max
    inset = ax.inset_axes([0.07, 0.5, 0.36, 0.42])
    inset.plot(np.linspace(10.0, t_max, 500), smooth_count(np.linspace(10.0, t_max, 500)),
               color=ORANGE, linewidth=2.4, solid_capstyle="butt")
    inset.step(np.concatenate(([10.0], zeros[inside], [t_max])),
               np.concatenate(([0], np.arange(1, inside.sum() + 1), [inside.sum()])),
               where="post", color=BLUE, linewidth=1.2)
    inset.set_xlim(10, t_max)
    inset.set_ylim(0, None)
    inset.tick_params(labelsize=7)
    inset.set_title("T up to 100", loc="left", fontsize=8, color=INK_SECONDARY)
    save(fig, "count_staircase")


def figure_spacings():
    t = read_csv(ZEROS_CSV)["t"]
    gap = np.diff(t)
    middle = 0.5 * (t[1:] + t[:-1])
    # Mean spacing near t is 2 pi / ln(t / 2 pi), so s has mean 1.
    s = gap * np.log(middle / (2.0 * np.pi)) / (2.0 * np.pi)

    # Summary numbers quoted in the README.
    smallest = np.argsort(s)[:2]
    print(f"spacings: {len(s):,} gaps, mean {s.mean():.4f}, std {s.std():.4f} "
          f"(Wigner surmise std {np.sqrt(3 * np.pi / 8 - 1):.4f})")
    for i in smallest:
        print(f"  small spacing {s[i]:.4f} at t = {middle[i]:.2f} (raw gap {gap[i]:.5f}, between zeros "
              f"{t[i]:.6f} and {t[i + 1]:.6f})")

    fig, ax = plt.subplots()
    ax.hist(s, bins=np.linspace(0, 3.5, 71), density=True, color=BLUE, alpha=0.55,
            edgecolor="white", linewidth=0.8, label=f"Zeros found")
    grid = np.linspace(0, 3.5, 400)
    ax.plot(grid, 32.0 / np.pi**2 * grid**2 * np.exp(-4.0 * grid**2 / np.pi), color=ORANGE, linewidth=1.8,
            label="Wigner surmise")
    ax.legend(loc="upper right")
    ax.set_xlabel("Normalised spacing s")
    ax.set_ylabel("Zero density")
    ax.set_title("Spacings between consecutive zeros (1 = average spacing)",
                 loc="left", fontsize=10, color=INK)
    ax.set_xlim(0, 3.5)
    save(fig, "spacings")


if __name__ == "__main__":
    if not os.path.exists(ZEROS_CSV):
        sys.exit(f"{ZEROS_CSV} not found: run `make run N=1000000` first")
    figure_z_curve()
    figure_term_errors()
    figure_count_staircase()
    figure_spacings()
