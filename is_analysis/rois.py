"""Export the segmentation as ImageJ ROIs, so every result can be looked at in Fiji.

One .zip per image, holding, for every z plane:

    z04_T7        outline of T cell 7 in plane 4
    z04_Tumour12  outline of tumour cell 12
    z04_IS01      the contact of synapse 1, i.e. what was measured as its area

Fiji reads the file with ROI Manager > More >> Open..., and "Show All with labels" draws the
names on the image. The numbers are the `T_Cell`, `Tumour_Cell` and `IS` columns of
results/synapse_table.csv.
"""
import numpy as np
import roifile
from skimage.measure import find_contours

from .segmentation import T_CELL, UNKNOWN

MIN_VERTICES = 4


def write(path, bodies, synapses, classes):
    """Write one ROI zip for an image. Returns the number of ROIs written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()  # roiwrite appends, and a stale file would keep old outlines
    rois = []
    for label, (kind, _) in sorted(classes.items()):
        tag = {T_CELL: "T", UNKNOWN: "Unclear"}.get(kind, "Tumour")
        rois += _outlines(bodies == label, f"{tag}{label}")
    for index in range(1, int(synapses.max()) + 1):
        rois += _outlines(synapses == index, f"IS{index:02d}")
    if rois:
        roifile.roiwrite(path, rois)
    return len(rois)


def _outlines(mask, name):
    """One polygon ROI per z plane the object appears in."""
    out = []
    for z, plane in enumerate(mask):
        if not plane.any():
            continue
        parts = [c for c in find_contours(np.pad(plane.astype(float), 1), 0.5)
                 if len(c) >= MIN_VERTICES]
        for part, contour in enumerate(parts):
            points = np.column_stack([contour[:, 1], contour[:, 0]]) - 1.0  # (x, y), undo the pad
            suffix = "" if len(parts) == 1 else chr(ord("a") + part)  # a cell can be cut in two
            roi = roifile.ImagejRoi.frompoints(
                points, name=f"z{z + 1:02d}_{name}{suffix}", z=z + 1)
            roi.stroke_color = _colour(name)
            out.append(roi)
    return out


def _colour(name):
    if name.startswith("IS"):
        return b"\xff\xff\xff\xff"  # white
    if name.startswith("T"):
        return b"\xff\x00\xe5\xff"  # cyan, as in the QC figures
    return b"\xff\xff\xd4\x00"      # yellow
