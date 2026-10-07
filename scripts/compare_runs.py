"""Put the automatic and the hand-placed segmentation of every image side by side.

    python scripts/compare_runs.py

Reads the two runs' QC projections and tables and writes one figure per image, plus a summary
of how many cells and synapses each run found. Useful for checking a curation pass, and the
record of what the automatic detection got wrong on this data.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def counts(table):
    """{image: (synapses, T cells, tumour cells)} from a synapse table."""
    if not table.exists():
        return {}
    df = pd.read_csv(table)
    return {image: (len(g), int(g.N_T_Cells.iloc[0]), int(g.N_Tumour_Cells.iloc[0]))
            for image, g in df.groupby("Image")}


def describe(count):
    if count is None:
        return "no cell found"
    synapses, t_cells, tumour = count
    return "{} T cells, {} tumour cells, {} synapse{}".format(
        t_cells, tumour, synapses, "" if synapses == 1 else "s")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--auto", default="results", help="output folder of the automatic run")
    ap.add_argument("--curated", default="results_curated", help="output folder of the hand run")
    ap.add_argument("--out", default="docs/auto_vs_curated")
    ap.add_argument("--dpi", type=int, default=95)
    a = ap.parse_args()

    auto, curated, out = Path(a.auto), Path(a.curated), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    auto_counts, curated_counts = counts(auto / "synapse_table.csv"), \
        counts(curated / "synapse_table.csv")

    rows = []
    for figure in sorted((curated / "segmentation_qc").glob("*.png")):
        name = figure.stem
        other = auto / "segmentation_qc" / figure.name
        if not other.exists():
            print(f"{name}: no automatic run, skipped")
            continue
        left, right = auto_counts.get(name), curated_counts.get(name)
        fig, axes = plt.subplots(1, 2, figsize=(14, 7.2))
        for ax, image, title in zip(axes, (other, figure),
                                    ("AUTOMATIC  -  " + describe(left),
                                     "HAND-PLACED  -  " + describe(right))):
            ax.imshow(mpimg.imread(image))
            ax.set_title(title, fontsize=12, pad=8)
            ax.axis("off")
        fig.text(0.5, 0.045, f"{name}.  cyan = T cell (GFP+), yellow = tumour cell, "
                             "white ring = a measured synapse.", ha="center", fontsize=10.5)
        fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.10, wspace=0.02)
        fig.savefig(out / f"{name}.png", dpi=a.dpi)
        plt.close(fig)
        rows.append(dict(Image=name,
                         Auto_Synapses=left[0] if left else 0,
                         Curated_Synapses=right[0] if right else 0,
                         Auto_T=left[1] if left else 0, Curated_T=right[1] if right else 0,
                         Auto_Tumour=left[2] if left else 0,
                         Curated_Tumour=right[2] if right else 0))

    summary = pd.DataFrame(rows)
    summary["Difference"] = summary.Curated_Synapses - summary.Auto_Synapses
    summary = summary.sort_values("Difference", ascending=False)
    summary.to_csv(out / "summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"\n{len(summary)} figures and summary.csv in {out}")


if __name__ == "__main__":
    main()
