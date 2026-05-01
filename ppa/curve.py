#!/usr/bin/env python3

import itertools
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":          10,
    "axes.linewidth":     0.8,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "xtick.major.size":   4,
    "ytick.major.size":   4,
    "legend.frameon":     True,
    "legend.framealpha":  0.9,
    "legend.edgecolor":   "#cccccc",
    "legend.fontsize":    8,
    "figure.dpi":         150,
})

# Data 
X = np.array([48.375, 0.376, 0.74425, 0.0754])
Y = np.array([0.625,  0.015, 0.025,   0.0045])

labels = [
    "Whole Plasma",
    "Whole CSF",
    "Plasma SEC Fr. 7–10",
    "CSF SEC Fr. 7–10",
]

# (x_offset_pts, y_offset_pts, ha)
annotation_offsets = {
    "Whole Plasma":        ( 10, -18, "left"),
    "Whole CSF":           (-10,  20, "right"),
    "Plasma SEC Fr. 7–10": ( 10, -18, "left"),
    "CSF SEC Fr. 7–10":    ( 10, -14, "left"),
}

grouped_points = {
    "Plasma 1": [(53.0,0.671215875),(12.02,0.21423457),(2.168,0.057326036),
                 (0.524,0.019215707),(0.167,0.007969149),(0.0947,0.005149721)],
    "Plasma 2": [(56.1,0.701235343),(8.65,0.166306177),(2.22,0.058381462),
                 (0.476,0.017846011),(0.149,0.007299423),(0.0828,0.00464403)],
    "Plasma 3": [(51,0.651634166),(22.4,0.345917857),(2.4,0.061992004),
                 (0.516,0.018989501),(0.0964,0.005220729),(0.0795,0.004500903)],
    "Plasma 4": [(51.6,0.65752695),(5.65,0.119822609),(2.184,0.057651397),
                 (0.488,0.018191304),(0.0895,0.004930662),(0.0753,0.004316743)],
    "Tagged EVs Only": [(0.0827, 0.004639713)],
}

# Fit 
def power_law(x, a, b):
    return a * x**b

params, _ = curve_fit(power_law, X, Y)
a, b = params
Y_pred = power_law(X, a, b)
r2 = r2_score(Y, Y_pred)

x_range = np.logspace(np.log10(min(X) * 0.4), np.log10(max(X) * 2.5), 300)
y_fit   = power_law(x_range, a, b)

# Colorblind-friendly palette (Wong 2011) 
WONG = {
    "blue":   "#0072B2",
    "orange": "#E69F00",
    "green":  "#009E73",
    "purple": "#CC79A7",
    "sky":    "#56B4E9",
    "black":  "#000000",
}

group_styles = {
    "Plasma 1":        dict(color=WONG["orange"], marker="o"),
    "Plasma 2":        dict(color=WONG["green"],  marker="s"),
    "Plasma 3":        dict(color=WONG["purple"], marker="^"),
    "Plasma 4":        dict(color=WONG["sky"],    marker="D"),
    "Tagged EVs Only": dict(color=WONG["black"],  marker="*"),
}

# Figure 
fig, ax = plt.subplots(figsize=(6.5, 5))

# Fit line
ax.plot(x_range, y_fit,
        color="#c0392b", linewidth=1.4, linestyle="--", zorder=2,
        label=f"Power-law fit: $Y = {a:.4f}\\,X^{{{b:.3f}}}$\n$R^2 = {r2:.4f}$")

# Reference points
ax.scatter(X, Y,
           color=WONG["blue"], marker="o", s=55, zorder=4,
           edgecolors="white", linewidths=0.5,
           label="Reference biofluid concentrations")

# Sample groups
for group_name, points in grouped_points.items():
    xv, yv = zip(*points)
    sty = group_styles[group_name]
    ax.scatter(xv, yv,
               label=group_name,
               color=sty["color"], marker=sty["marker"],
               s=40, alpha=0.85, zorder=3,
               edgecolors="white", linewidths=0.4)
    
# Annotations 
for x, y, lbl in zip(X, Y, labels):
    ox, oy, ha = annotation_offsets[lbl]
    ax.annotate(
        lbl, xy=(x, y),
        xytext=(ox, oy),
        textcoords="offset points",
        fontsize=7.5, color="#2c2c2c",
        ha=ha, va="center",
        arrowprops=dict(
            arrowstyle="-",
            color="#888888",
            lw=0.6,
            shrinkA=4,
            shrinkB=3,
        ),
    )

# Axes 
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(left=0.04, right=130)
ax.set_ylim(bottom=0.002, top=1.5)

ax.set_xlabel("Protein Concentration (mg mL$^{-1}$)", fontsize=10, labelpad=6)
ax.set_ylabel("Proteinase K Concentration (mg mL$^{-1}$)", fontsize=10, labelpad=6)

def log_formatter(val, pos):
    exp = int(round(np.log10(val)))
    if abs(val - 10**exp) < 1e-9 * val:
        return f"$10^{{{exp}}}$"
    return ""

for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_formatter(ticker.FuncFormatter(log_formatter))
    axis.set_minor_locator(ticker.LogLocator(subs=np.arange(2, 10)))
    axis.set_minor_formatter(ticker.NullFormatter())

ax.tick_params(axis="both", which="minor", length=2.5, width=0.6)

# Legend
leg = ax.legend(
    loc="upper left",
    borderpad=0.7,
    labelspacing=0.45,
    handlelength=1.6,
    handleheight=0.9,
    handletextpad=0.5,
)
leg.get_frame().set_linewidth(0.6)

# Save 
plt.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.12)
plt.savefig("protein_conc_pk_curve.png", dpi=300, facecolor="white")