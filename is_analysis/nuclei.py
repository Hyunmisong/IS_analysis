"""Nucleus detection from the DAPI channel.

Two interchangeable detectors, both returning a 3D label image with one label per nucleus:

- `cellpose_nuclei` (default) runs the Cellpose "nuclei" model, plane by plane, and stitches the
  planes into 3D objects. It needs `pip install cellpose`.
- `watershed_nuclei` needs nothing beyond scikit-image: threshold, then split touching nuclei on
  the distance map. On low-photon images its mask breaks into pieces, which costs a nucleus its
  seed or gives it several, so the cell it belongs to is merged into its neighbour or cut in two.
  It is kept as the dependency-free fallback.


Only the *seeds* matter downstream: the cell body is built from all three channels in
`segmentation.cell_bodies`, and the nucleus label is used to place and classify the cell. A
nucleus mask that is a little tight is therefore harmless, a missing or duplicated one is not.
"""
import math

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed

from .segmentation import ball, smooth

DIAMETER_UM = 8.0
SMOOTH_UM = 1.0
STITCH = 0.25
MIN_FRACTION = 0.075
_MODEL = None


def smallest_kept_um3(diameter_um, fraction):
    """Volume below which a detected object is thrown away, as a fraction of the volume of a
    sphere of `diameter_um`, so the cut-off follows the cell type rather than being absolute.

    It defaults to 0.075 - 20 um3 for an 8 um nucleus, a ball of 3.4 um across - which is a floor
    against single specks of noise and nothing more. A real cut-off does more harm than good on
    the 2026-06-04 images: a fifth of the detections fall below 0.2 of the expected volume,
    but checking them one by one they are not debris: some are nuclei Cellpose simply drew a
    little small, and the rest are cells sitting above or below another cell - their outlines
    overlap in x and y but they are five to nine planes apart in z, so they are not pieces of
    one nucleus either. Removing them cost 7 of 60 synapses and *raised* the largest contact
    from 194 to 263 um2, because a cell whose only seed is gone has no marker left and the
    watershed hands its body to the neighbour, turning their shared border into a false contact.

    Raise it only for data where small debris really is present, and check the synapse count
    before and after.
    """
    return fraction * math.pi / 6 * diameter_um ** 3


def segment(dapi, voxel, method="cellpose", min_fraction=MIN_FRACTION, min_sep_um=4.0,
            diameter_um=DIAMETER_UM, smooth_um=SMOOTH_UM):
    if method == "cellpose":
        labels = cellpose_nuclei(dapi, voxel, diameter_um, smooth_um)
    else:
        labels = watershed_nuclei(dapi, voxel, min_sep_um=min_sep_um)
    return _renumber(labels, voxel, smallest_kept_um3(diameter_um, min_fraction))


def cellpose_nuclei(dapi, voxel, diameter_um=DIAMETER_UM, smooth_um=SMOOTH_UM, stitch=STITCH):
    """Cellpose on a smoothed DAPI stack, 2D per plane, stitched in z.

    The raw channel is far too grainy for the model - on these images it finds almost nothing
    until the shot noise is smoothed away, so `smooth_um` is applied first. Plane-by-plane with
    stitching beats Cellpose's own 3D mode here: the stacks are strongly anisotropic (1 um in z
    against 0.2 um in x and y) and the 3D mode returns hollow, ring-shaped objects, five times
    slower.
    """
    model = _cellpose_model()
    labels, _, _ = model.eval(smooth(dapi, voxel, smooth_um), channels=[0, 0], z_axis=0,
                              diameter=diameter_um / voxel[2], stitch_threshold=stitch,
                              normalize=True)
    return _renumber(np.asarray(labels, np.int32), voxel)


def watershed_nuclei(dapi, voxel, min_sep_um=4.0, sigma_um=0.4):
    """Threshold the DAPI channel and split touching nuclei on the distance map."""
    sm = smooth(dapi, voxel, sigma_um)
    mask = ndi.binary_fill_holes(ndi.binary_closing(sm > threshold_otsu(sm), ball(0.5, voxel)))
    dist = ndi.distance_transform_edt(mask, sampling=voxel)
    seeds = ndi.label(_peaks(smooth(dist, voxel, min_sep_um / 4), mask, voxel, min_sep_um))[0]
    return _renumber(watershed(-dist, seeds, mask=mask), voxel)


def _cellpose_model():
    global _MODEL
    if _MODEL is None:
        try:
            from cellpose import models
        except ImportError:
            raise SystemExit("--nuclei cellpose needs the cellpose package: pip install cellpose\n"
                             "(or run with --nuclei watershed, which needs nothing extra)")
        _MODEL = models.CellposeModel(gpu=False, model_type="nuclei")
    return _MODEL


def _peaks(dist, mask, voxel, min_sep_um):
    """One seed per maximum of the distance map, maxima closer than min_sep_um merged.

    The map is smoothed by min_sep_um / 4 first: a large nucleus has a lumpy map and would
    otherwise raise several maxima and be cut into pieces.
    """
    return mask & (dist >= ndi.maximum_filter(dist, footprint=ball(min_sep_um / 2, voxel))) \
        & (dist > 0)


def _renumber(labels, voxel, min_volume_um3=0.0):
    """Drop anything below min_volume_um3 and make the remaining labels 1..n."""
    keep = np.bincount(labels.ravel()) * np.prod(voxel) > min_volume_um3
    keep[0] = False
    out = np.zeros(labels.shape, np.int32)
    for new, old in enumerate(np.flatnonzero(keep), start=1):
        out[labels == old] = new
    return out
