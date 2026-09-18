#!/usr/bin/env python3
"""Reproduce the TCGA-LUAD PLVAP/PECAM1 survival analysis from UCSC Xena.

The script downloads the phenotype, overall-survival, PLVAP, and PECAM1 data;
constructs log2((PLVAP TPM + 0.001)/(PECAM1 TPM + 0.001)); divides patients at
the median; runs Kaplan-Meier, log-rank, and Cox analyses; and exports the
patient-level data, statistics, PNG, SVG, PDF, and provenance metadata.

Example:
    python reproduce_LUAD_PLVAP_PECAM1_ratio_from_xena.py
"""

from argparse import ArgumentParser
from datetime import datetime, timezone
from pathlib import Path
import json
import math
import platform
import warnings

import lifelines
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xenaPython as xena


HUB = "https://toil.xenahubs.net"
PHENOTYPE_DATASET = "TcgaTargetGTEX_phenotype.txt"
SURVIVAL_DATASET = "TCGA_survival_data"
EXPRESSION_DATASET = "tcga_RSEM_gene_tpm"
GENES = ["PLVAP", "PECAM1"]
DAYS_PER_MONTH = 365.25 / 12
FONT_SCALE = 2.0
COLORS = {"Low": "#0571B0", "High": "#CA0020"}
DEFAULT_PREFIX = "LUAD_PLVAP_PECAM1_ratio"

plt.rcParams["svg.fonttype"] = "none"


def fs(size):
    return size * FONT_SCALE


def decode_category(field, values):
    codes = xena.field_codes(HUB, PHENOTYPE_DATASET, [field])[0]["code"].split(
        "\t"
    )
    return [None if str(value) == "NaN" else codes[int(float(value))] for value in values]


def fetch_data():
    """Download and merge the LUAD phenotype, survival, and expression data."""
    phenotype_samples = xena.dataset_samples(HUB, PHENOTYPE_DATASET, None)
    phenotype_fields = ["_study", "_sample_type", "primary disease or tissue"]
    _, phenotype_values = xena.dataset_probe_values(
        HUB, PHENOTYPE_DATASET, phenotype_samples, phenotype_fields
    )
    decoded = {
        field: decode_category(field, values)
        for field, values in zip(phenotype_fields, phenotype_values)
    }
    luad_primary = [
        sample
        for index, sample in enumerate(phenotype_samples)
        if decoded["_study"][index] == "TCGA"
        and decoded["_sample_type"][index] == "Primary Tumor"
        and decoded["primary disease or tissue"][index] == "Lung Adenocarcinoma"
    ]

    survival_samples = xena.dataset_samples(HUB, SURVIVAL_DATASET, None)
    _, survival_values = xena.dataset_probe_values(
        HUB, SURVIVAL_DATASET, survival_samples, ["OS", "OS.time"]
    )
    survival_map = {
        sample: (event, days)
        for sample, event, days in zip(survival_samples, *survival_values)
    }

    expression_results = xena.dataset_gene_probe_avg(
        HUB, EXPRESSION_DATASET, luad_primary, GENES
    )
    expression = {
        result["gene"]: np.asarray(result["scores"][0], dtype=float)
        for result in expression_results
    }

    rows = []
    for index, sample in enumerate(luad_primary):
        event_raw, days_raw = survival_map.get(sample, ("NaN", "NaN"))
        plvap = float(expression["PLVAP"][index])
        pecam1 = float(expression["PECAM1"][index])
        survival_present = str(event_raw) != "NaN" and str(days_raw) != "NaN"
        expression_present = np.isfinite(plvap) and np.isfinite(pecam1)
        included = survival_present and expression_present

        if included:
            event = int(float(event_raw))
            days = float(days_raw)
            if event not in (0, 1) or days < 0:
                included = False
        else:
            event = np.nan
            days = np.nan

        reasons = []
        if not survival_present:
            reasons.append("Missing OS event or time")
        if not expression_present:
            reasons.append("Missing PLVAP or PECAM1 expression")
        if survival_present and (event not in (0, 1) or days < 0):
            reasons.append("Invalid OS event or time")

        rows.append(
            {
                "sample_id": sample,
                "cancer": "LUAD",
                "sample_type": "Primary Tumor",
                "PLVAP_log2_tpm_plus_0.001": plvap,
                "PECAM1_log2_tpm_plus_0.001": pecam1,
                "PLVAP_PECAM1_log2_ratio": plvap - pecam1,
                "OS_event": event if included else np.nan,
                "OS_days": days if included else np.nan,
                "OS_months": days / DAYS_PER_MONTH if included else np.nan,
                "analysis_included": included,
                "exclusion_reason": "; ".join(reasons),
            }
        )
    return pd.DataFrame(rows)


def analyze(data):
    """Apply the median split and calculate survival statistics."""
    analyzed = data.loc[data["analysis_included"]].copy()
    score = analyzed["PLVAP_PECAM1_log2_ratio"]
    cutoff = float(score.median())
    analyzed["PLVAP_PECAM1_ratio_group"] = np.where(
        score < cutoff, "Low", "High"
    )

    low = analyzed["PLVAP_PECAM1_ratio_group"] == "Low"
    high = analyzed["PLVAP_PECAM1_ratio_group"] == "High"
    if low.sum() != high.sum():
        raise ValueError(
            "The median falls on tied ratio values, so an exact 50/50 split is "
            "not possible with a simple threshold. Specify a tie rule before proceeding."
        )

    logrank = logrank_test(
        analyzed.loc[low, "OS_months"],
        analyzed.loc[high, "OS_months"],
        event_observed_A=analyzed.loc[low, "OS_event"],
        event_observed_B=analyzed.loc[high, "OS_event"],
    )
    cox_data = analyzed[["OS_months", "OS_event"]].copy()
    cox_data["high_group"] = high.astype(int).to_numpy()
    cox = CoxPHFitter().fit(
        cox_data, duration_col="OS_months", event_col="OS_event"
    )

    km = {}
    median_os = {}
    for group in ["Low", "High"]:
        mask = analyzed["PLVAP_PECAM1_ratio_group"] == group
        km[group] = KaplanMeierFitter(label=group).fit(
            analyzed.loc[mask, "OS_months"],
            analyzed.loc[mask, "OS_event"],
        )
        median_os[group] = float(km[group].median_survival_time_)

    ci = np.exp(cox.confidence_intervals_.loc["high_group"]).to_numpy()
    stats = {
        "analysis": "TCGA-LUAD overall survival by PLVAP/PECAM1 expression ratio",
        "ratio_definition": "log2((PLVAP TPM + 0.001)/(PECAM1 TPM + 0.001))",
        "grouping": "median ratio cutoff; score below median is Low",
        "PLVAP_PECAM1_log2_ratio_median_cutoff": cutoff,
        "n_total_luad_primary": int(len(data)),
        "n_analyzed": int(len(analyzed)),
        "n_low": int(low.sum()),
        "n_high": int(high.sum()),
        "events": int(analyzed["OS_event"].sum()),
        "median_os_low_months": median_os["Low"],
        "median_os_high_months": median_os["High"],
        "logrank_chi_square": float(logrank.test_statistic),
        "logrank_p": float(logrank.p_value),
        "hazard_ratio_high_vs_low": float(np.exp(cox.params_["high_group"])),
        "hazard_ratio_ci95_low": float(ci[0]),
        "hazard_ratio_ci95_high": float(ci[1]),
        "hazard_ratio_p": float(cox.summary.loc["high_group", "p"]),
        "days_per_month": DAYS_PER_MONTH,
    }

    result = data.copy()
    result["PLVAP_PECAM1_ratio_group"] = ""
    result.loc[analyzed.index, "PLVAP_PECAM1_ratio_group"] = analyzed[
        "PLVAP_PECAM1_ratio_group"
    ]
    return result, analyzed, km, stats


def make_plot(analyzed, km, stats, output_dir, prefix):
    max_month = float(analyzed["OS_months"].max())
    tick_step = 50 if max_month > 160 else 25
    max_tick = int(math.ceil(max_month / tick_step) * tick_step)
    ticks = np.arange(0, max_tick + 0.1, tick_step)

    fig = plt.figure(figsize=(10.0, 10.0), facecolor="white")
    ax = fig.add_axes([0.20, 0.23, 0.75, 0.47])
    for group in ["Low", "High"]:
        km[group].plot_survival_function(
            ax=ax,
            ci_show=True,
            color=COLORS[group],
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
        "TCGA-LUAD survival by\nPLVAP/PECAM1 expression ratio",
        ha="left",
        va="top",
        fontsize=fs(15),
        weight="bold",
        linespacing=0.92,
    )
    fig.text(
        0.10,
        0.815,
        "log2 expression ratio;\nmedian cutoff (50/50)",
        ha="left",
        va="top",
        fontsize=fs(9.5),
        color="#4B5563",
    )
    ax.text(
        0.96,
        0.94,
        f"Log-rank p = {stats['logrank_p']:.4f}\n"
        f"HR = {stats['hazard_ratio_high_vs_low']:.2f} "
        f"(95% CI {stats['hazard_ratio_ci95_low']:.2f}-"
        f"{stats['hazard_ratio_ci95_high']:.2f})",
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
        ["Low PLVAP/PECAM1\nratio", "High PLVAP/PECAM1\nratio"],
        title=None,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.18),
        ncol=2,
        columnspacing=2.0,
        handlelength=2.0,
        fontsize=fs(10),
    )
    ax.tick_params(axis="both", labelsize=fs(10))
    ax.set_xticks(ticks)
    fig.text(
        0.10,
        0.018,
        "Data: UCSC Xena Toil RNA-seq recompute and TCGA\n"
        "Pan-Cancer survival endpoints.",
        fontsize=fs(7.8),
        color="#6B7280",
        ha="left",
        va="bottom",
    )

    plot_stem = output_dir / f"{prefix}_KM_2x_clean_taller_yaxis"
    for extension in ["png", "svg", "pdf"]:
        fig.savefig(plot_stem.with_suffix(f".{extension}"), dpi=300)
    plt.close(fig)


def write_outputs(result, stats, output_dir, prefix):
    result.to_csv(
        output_dir / f"{prefix}_patient_data.csv",
        index=False,
        float_format="%.8g",
    )
    with (output_dir / f"{prefix}_statistics.json").open("w") as handle:
        json.dump(stats, handle, indent=2)
    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "hub": HUB,
        "phenotype_dataset": PHENOTYPE_DATASET,
        "survival_dataset": SURVIVAL_DATASET,
        "expression_dataset": EXPRESSION_DATASET,
        "genes": GENES,
        "sample_filter": {
            "study": "TCGA",
            "sample_type": "Primary Tumor",
            "primary_disease_or_tissue": "Lung Adenocarcinoma",
        },
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "lifelines": lifelines.__version__,
    }
    with (output_dir / f"{prefix}_provenance.json").open("w") as handle:
        json.dump(provenance, handle, indent=2)


def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output directory (default: directory containing this script)",
    )
    parser.add_argument(
        "--prefix",
        default=DEFAULT_PREFIX,
        help=f"Output filename prefix (default: {DEFAULT_PREFIX})",
    )
    return parser.parse_args()


def main():
    warnings.filterwarnings("ignore", category=FutureWarning)
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print("Downloading phenotype, survival, PLVAP, and PECAM1 data...", flush=True)
    source = fetch_data()
    print("Analyzing survival and generating outputs...", flush=True)
    result, analyzed, km, stats = analyze(source)
    write_outputs(result, stats, args.output_dir, args.prefix)
    make_plot(analyzed, km, stats, args.output_dir, args.prefix)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
