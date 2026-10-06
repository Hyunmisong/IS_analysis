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

from .segmentation import T_CELL, TUMOUR, ball, smooth

SEED_RADIUS_UM = 1.5
KINDS = {"t": T_CELL, "tcell": T_CELL, "t_cell": T_CELL,
         "tumour": TUMOUR, "tumor": TUMOUR}


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
