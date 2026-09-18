#!/usr/bin/env python3
"""Regenerate the clean LUAD plot with a 25% taller y-axis and clear legend."""

from pathlib import Path
import json
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter


HERE = Path(__file__).resolve().parent
FONT_SCALE = 2.0

# Keep SVG text as text so the scaled font sizes remain directly inspectable.
plt.rcParams["svg.fonttype"] = "none"


def fs(original_size):
    """Return an original point size scaled by the requested factor."""
    return original_size * FONT_SCALE


def at_risk(frame, times):
    return [int((frame["OS_months"] >= time).sum()) for time in times]


data = pd.read_csv(HERE / "LUAD_PLVAP_POSTN_patient_data.csv")
analyzed = data.loc[data["analysis_included"] == True].copy()  # noqa: E712
with (HERE / "LUAD_PLVAP_POSTN_statistics.json").open() as handle:
    stats = json.load(handle)

km = {}
for group in ["Low", "High"]:
    group_data = analyzed.loc[analyzed["group"] == group]
    km[group] = KaplanMeierFitter(
        label=f"{group} PLVAP + POSTN\nsignature"
    ).fit(
        group_data["OS_months"], group_data["OS_event"]
    )

colors = {"Low": "#0571B0", "High": "#CA0020"}
max_month = float(analyzed["OS_months"].max())
tick_step = 50 if max_month > 160 else 25
max_tick = int(math.ceil(max_month / tick_step) * tick_step)
ticks = np.arange(0, max_tick + 0.1, tick_step)

# The number-at-risk panel is intentionally omitted in this clean version.
# Increasing the canvas height from 8 to 10 inches makes the plotting axis
# exactly 25% taller while preserving its normalized layout and x dimension.
fig = plt.figure(figsize=(10.0, 10.0), facecolor="white")
ax = fig.add_axes([0.20, 0.23, 0.75, 0.47])

for group in ["Low", "High"]:
    km[group].plot_survival_function(
        ax=ax,
        ci_show=True,
        color=colors[group],
        linewidth=2.2,
        ci_alpha=0.13,
        show_censors=True,
        censor_styles={"ms": 4, "marker": "+", "mew": 1},
    )

ax.set_xlim(0, max_tick)
ax.set_ylim(0, 1.02)
ax.set_ylabel("Overall survival\nprobability", fontsize=fs(11), linespacing=0.90)
ax.set_xlabel("Time (months)", fontsize=fs(11), labelpad=fs(8))
fig.text(
    0.10,
    0.965,
    "TCGA-LUAD survival by\nPLVAP + POSTN signature",
    ha="left",
    va="top",
    fontsize=fs(15),
    weight="bold",
    linespacing=0.92,
)
fig.text(
    0.10,
    0.815,
    "Unscaled PC1 signature;\nmedian cutoff (50/50)",
    ha="left",
    va="top",
    fontsize=fs(9.5),
    color="#4B5563",
)
annotation = (
    f"Log-rank p = {stats['logrank_p']:.4f}\n"
    f"HR = {stats['hazard_ratio_high_vs_low']:.2f} "
    f"(95% CI {stats['hazard_ratio_ci95_low']:.2f}-{stats['hazard_ratio_ci95_high']:.2f})"
)
ax.text(
    0.96,
    0.94,
    annotation,
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=fs(9.5),
)
ax.grid(False)
ax.spines[["top", "right"]].set_visible(False)
legend_handles, _ = ax.get_legend_handles_labels()
ax.legend(
    legend_handles,
    ["Low PLVAP + POSTN\nsignature", "High PLVAP + POSTN\nsignature"],
    title=None,
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.18),
    ncol=2,
    columnspacing=1.4,
    handlelength=2.0,
    fontsize=fs(10),
)
ax.tick_params(axis="both", labelsize=fs(10))
ax.set_xticks(ticks)

fig.text(
    0.10,
    0.018,
    "Data: UCSC Xena Toil RNA-seq recompute and TCGA\nPan-Cancer survival endpoints.",
    fontsize=fs(7.8),
    color="#6B7280",
    ha="left",
    va="bottom",
)

for extension in ["png", "svg", "pdf"]:
    fig.savefig(
        HERE / f"LUAD_PLVAP_POSTN_KM_2x_clean_signature_labels_taller_yaxis.{extension}",
        dpi=300,
    )
plt.close(fig)
