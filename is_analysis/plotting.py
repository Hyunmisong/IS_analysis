"""QC figures (per image, per synapse) and the final comparison plot."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import find_boundaries

from .segmentation import T_CELL

COLOURS = {T_CELL: "#00e5ff", "Tumour": "#ffd400"}


def rgb(stack, gamma=0.6):
    """DAPI blue, GFP green, mCherry red maximum-intensity projection."""
    dapi, gfp, mcherry = (_norm(c.max(axis=0)) ** gamma for c in stack)
    return np.dstack([mcherry, gfp, dapi])


def _norm(img, low=1, high=99.8):
    lo, hi = np.percentile(img, [low, high])
    return np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)


def plot_image_qc(path, stack, bodies, classes, rows, voxel, title):
    """Projection with the cell outlines coloured by cell type and every synapse marked."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(rgb(stack))
    flat = {label: (bodies == label).any(axis=0) for label in classes}  # a small cell can be
    for label, (kind, _) in classes.items():                             # hidden by a larger one
        ax.contour(find_boundaries(flat[label], mode="inner"), levels=[0.5],
                   colors=COLOURS[kind], linewidths=0.8)
        y, x = np.argwhere(flat[label]).mean(axis=0)
        ax.text(x, y, str(label), color=COLOURS[kind], fontsize=7, ha="center", va="center")
    for r in rows:
        y, x = _centre(flat[r["T_Cell"]], flat[r["Tumour_Cell"]])
        ax.plot(x, y, "o", mfc="none", mec="w", ms=14, mew=1.2)
        ax.text(x + 9, y, f"{r['Contact_Area_um2']:.1f} um2", color="w", fontsize=7, va="center")
    _scalebar(ax, voxel[2], bodies.shape[1:])
    ax.set_title(title, fontsize=9)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_synapse(path, stack, bodies, row, voxel, title):
    """The contact plane: merge / mCherry / cell bodies, cropped around the cell pair."""
    path.parent.mkdir(parents=True, exist_ok=True)
    t, tumour = row["T_Cell"], row["Tumour_Cell"]
    pair = (bodies == t) | (bodies == tumour)
    z = int(round(row["Z_Contact"]))
    yx = np.argwhere(pair.any(axis=0))
    (y0, x0), (y1, x1) = yx.min(0) - 5, yx.max(0) + 6
    sl = (slice(max(y0, 0), y1), slice(max(x0, 0), x1))
    plane = tuple(c[z][sl] for c in stack)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9))
    axes[0].imshow(np.dstack([_norm(plane[2]), _norm(plane[1]), _norm(plane[0])]))
    axes[0].set_title("merge (R mCherry / G GFP / B DAPI)", fontsize=8)
    axes[1].imshow(plane[2], cmap="inferno")
    axes[1].set_title(f"mCherry  IS/surface = {row['mCherry_Enrichment']:.2f}", fontsize=8)
    axes[2].imshow((bodies[z][sl] == t) * 1 + (bodies[z][sl] == tumour) * 2, cmap="viridis",
                   vmin=0, vmax=2)
    axes[2].set_title(f"T (mid) vs tumour (bright), z={z}", fontsize=8)
    for ax in axes:
        ax.contour((bodies[z][sl] == t) & _dilate_plane(bodies[z][sl] == tumour), levels=[0.5],
                   colors="w", linewidths=1.3)
        _scalebar(ax, voxel[2], plane[0].shape)
        ax.axis("off")
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _dilate_plane(mask):
    from scipy import ndimage as ndi

    return ndi.binary_dilation(mask, np.ones((3, 3), bool))


def _centre(a, b):
    """Where to put the marker: the middle of the a/b border in projection."""
    border = a & _dilate_plane(b)
    return np.argwhere(border if border.any() else a).mean(axis=0)


def _scalebar(ax, dx_um, shape, length_um=5.0):
    n = length_um / dx_um
    y, x = shape[0] * 0.95, shape[1] * 0.05
    ax.plot([x, x + n], [y, y], "w-", lw=2.5)
    ax.text(x, y - shape[0] * 0.025, f"{length_um:g} um", color="w", fontsize=7)


def plot_groups(path, df, value, control, group_p, title):
    """Box plot of one metric per condition with the significance vs the control group."""
    groups = list(df.groupby("Condition"))
    fig, ax = plt.subplots(figsize=(1.5 * len(groups) + 2, 4.2))
    data = [g[value].dropna().values for _, g in groups]
    ax.boxplot(data, tick_labels=[c for c, _ in groups], showfliers=False)
    for i, values in enumerate(data, start=1):
        ax.plot(np.random.normal(i, 0.055, len(values)), values, "o", ms=3.5, alpha=0.6,
                color="#2b6cb0")
    top = max((v.max() for v in data if len(v)), default=1.0)
    for i, (condition, _) in enumerate(groups, start=1):
        p = group_p.get(condition, float("nan"))
        if condition != control and p == p:
            ax.text(i, top * 1.04, _stars(p), ha="center", fontsize=10)
    ax.set_ylabel(value.replace("_", " "))
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
