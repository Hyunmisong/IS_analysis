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
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed

from .segmentation import ball, smooth

DIAMETER_UM = 8.0
SMOOTH_UM = 1.0
STITCH = 0.25
_MODEL = None


def segment(dapi, voxel, method="cellpose", min_volume_um3=20.0, min_sep_um=4.0,
            diameter_um=DIAMETER_UM, smooth_um=SMOOTH_UM):
    if method == "cellpose":
        return cellpose_nuclei(dapi, voxel, diameter_um, smooth_um, min_volume_um3)
    return watershed_nuclei(dapi, voxel, min_volume_um3, min_sep_um)


def cellpose_nuclei(dapi, voxel, diameter_um=DIAMETER_UM, smooth_um=SMOOTH_UM,
                    min_volume_um3=20.0, stitch=STITCH):
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
    return _renumber(np.asarray(labels, np.int32), voxel, min_volume_um3)


def watershed_nuclei(dapi, voxel, min_volume_um3=20.0, min_sep_um=4.0, sigma_um=0.4):
    """Threshold the DAPI channel and split touching nuclei on the distance map."""
    sm = smooth(dapi, voxel, sigma_um)
    mask = ndi.binary_fill_holes(ndi.binary_closing(sm > threshold_otsu(sm), ball(0.5, voxel)))
    mask = remove_small_objects(mask, int(min_volume_um3 / np.prod(voxel)))
    dist = ndi.distance_transform_edt(mask, sampling=voxel)
    seeds = ndi.label(_peaks(smooth(dist, voxel, min_sep_um / 4), mask, voxel, min_sep_um))[0]
    return _renumber(watershed(-dist, seeds, mask=mask), voxel, min_volume_um3)


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


def _renumber(labels, voxel, min_volume_um3):
    """Drop everything below min_volume_um3 and make the remaining labels 1..n."""
    keep = np.bincount(labels.ravel()) * np.prod(voxel) >= min_volume_um3
    keep[0] = False
    out = np.zeros(labels.shape, np.int32)
    for new, old in enumerate(np.flatnonzero(keep), start=1):
        out[labels == old] = new
    return out
