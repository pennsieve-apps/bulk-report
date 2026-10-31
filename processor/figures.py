"""The report's own figure: per-library QC bars (PNG, 150 dpi)."""
from __future__ import annotations

from pathlib import Path

DPI = 150
BAR, FLAG, TEXT, TEXT2, GRID = "#2a78d6", "#e34948", "#0b0b0b", "#52514e", "#e4e3df"


def qc_bars(qc: dict, path: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2})
    rows = qc["rows"]
    ids = [str(r["sample_id"]) for r in rows]
    colours = [FLAG if r["flags"] else BAR for r in rows]
    if qc["source"] == "quant":
        panels = [("reads after trimming (M)", [(r["reads_after_trim"] or r["reads_in"]) / 1e6 for r in rows]),
                  ("mapping rate (%)", [100 * (r["mapping_rate"] or 0) for r in rows])]
    else:
        panels = [("counts per library (M)", [r["library_size"] / 1e6 for r in rows]),
                  ("genes detected (k)", [r["n_genes_detected"] / 1e3 for r in rows])]
    fig, axes = plt.subplots(1, 2, figsize=(max(6.0, 0.32 * len(ids) * 2 + 2), 3.0))
    for ax, (label, vals) in zip(axes, panels):
        ax.bar(range(len(ids)), vals, color=colours, width=0.7)
        ax.set_xticks(range(len(ids)), ids, rotation=60, ha="right", fontsize=7)
        ax.set_ylabel(label, color=TEXT)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path
