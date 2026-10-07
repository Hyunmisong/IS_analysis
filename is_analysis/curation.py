"""Hand-placed cell markers: one point per cell, plus what kind of cell it is.

Every failure the automatic detection makes on this data is a *seed* failure - one nucleus
detected as two, two nuclei detected as one, or the wrong cell type - and once two nuclei have
been merged, no later step gets them apart again (`docs/parameter-tuning.md`). The cell bodies
grown from correct seeds are good: on synthetic pairs the contact area comes out within 3-40 %
of the truth. So the part worth doing by hand is the seeds, and only the seeds.

A curation file is a three-column CSV in image pixel coordinates:

    x,y,type
    241,118,T
    263,131,Tumour

z is not needed. Each point is placed at the plane where the DAPI column under it is brightest,
which is where that nucleus sits, so points can all be clicked on one convenient plane.
"""
import csv
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from .segmentation import T_CELL, TUMOUR, UNKNOWN, background, ball, smooth

SEED_RADIUS_UM = 1.5
KINDS = {"t": T_CELL, "tcell": T_CELL, "t_cell": T_CELL,
         "tumour": TUMOUR, "tumor": TUMOUR,
         "unclear": UNKNOWN, "unknown": UNKNOWN, "?": UNKNOWN}
SNR_BOX_UM = 3.0


def file_for(directory, image_name):
    """The curation file for one image, or None if it has not been curated."""
    if not directory:
        return None
    path = Path(directory) / f"{image_name}.csv"
    return path if path.exists() else None


def read_points(path):
    """[(x_pixel, y_pixel, kind)] from a curation CSV."""
    points = []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            kind = KINDS.get(row["type"].strip().lower().replace(" ", ""))
            if kind is None:
                raise SystemExit(f"{path}: cell type must be T or Tumour, not {row['type']!r}")
            points.append((float(row["x"]), float(row["y"]), kind))
    if not points:
        raise SystemExit(f"{path}: no points")
    return points


def seeds(points, dapi, voxel, radius_um=SEED_RADIUS_UM):
    """(label image, {label: (kind, nan)}) from hand-placed points.

    Each point becomes a small ball at the brightest plane of its own DAPI column. The ball is
    only a marker for the watershed that grows the cell bodies, so its size does not matter as
    long as two neighbouring markers stay apart - which is why a point already taken by an
    earlier marker is never overwritten.
    """
    profile = smooth(dapi, voxel, 1.0)
    shape = dapi.shape
    marker = ball(radius_um, voxel)
    offsets = np.argwhere(marker) - (np.array(marker.shape) - 1) // 2
    labels = np.zeros(shape, np.int32)
    classes = {}
    for index, (x, y, kind) in enumerate(points, start=1):
        col, row = int(round(x)), int(round(y))
        if not (0 <= row < shape[1] and 0 <= col < shape[2]):
            raise SystemExit(f"point ({x}, {y}) is outside the {shape[2]} x {shape[1]} image")
        plane = int(np.argmax(profile[:, row, col]))
        here = offsets + np.array([plane, row, col])
        inside = np.all((here >= 0) & (here < shape), axis=1)
        here = here[inside]
        free = labels[tuple(here.T)] == 0
        labels[tuple(here[free].T)] = index
        classes[index] = (kind, float("nan"))
    return labels, classes


def write_points(path, nuclei, classes):
    """Write the automatic detection as a curation file, to be edited by hand."""
    path.parent.mkdir(parents=True, exist_ok=True)
    centres = ndi.center_of_mass(nuclei > 0, nuclei, sorted(classes))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "y", "type"])
        for (_, row, col), label in zip(centres, sorted(classes)):
            if np.isnan(row):
                continue
            writer.writerow([int(round(col)), int(round(row)), classes[label][0]])


def gfp_snr(points, dapi, gfp, voxel, box_um=SNR_BOX_UM):
    """GFP brightness of each hand-placed cell, in robust SDs above the image background.

    Measured in a `box_um` box at the brightest plane of the cell's own DAPI column - the same
    place the seed goes - so it answers "how green is this cell" without depending on a body
    having been grown for it yet.
    """
    smoothed = smooth(gfp, voxel, 0.5)
    bg, sd = background(gfp)
    profile = smooth(dapi, voxel, 1.0)
    half = [max(1, int(round(box_um / v))) for v in voxel]
    out = []
    for x, y, _ in points:
        col, row = int(round(x)), int(round(y))
        plane = int(np.argmax(profile[:, row, col]))
        box = smoothed[max(plane - half[0], 0):plane + half[0] + 1,
                       max(row - half[1], 0):row + half[1] + 1,
                       max(col - half[2], 0):col + half[2] + 1]
        out.append((float(box.mean()) - bg) / sd if sd else float("nan"))
    return out


def mark_unclear(points, snrs, low, high):
    """Set the type to Unclear where the GFP brightness falls in the [low, high) band.

    In that band the hand calls are about half T cell and half tumour cell at the same
    brightness, so the image is not deciding it. Such a cell keeps its point - without a seed its
    body would be absorbed by a neighbour and spoil *that* cell's boundary too - but it never
    forms a synapse.
    """
    marked = []
    for (x, y, kind), snr in zip(points, snrs):
        if low <= snr < high:
            kind = UNKNOWN
        marked.append((x, y, kind))
    return marked
