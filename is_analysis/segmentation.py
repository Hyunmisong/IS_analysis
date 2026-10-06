"""3D nuclei segmentation, T cell / tumour cell classification and cell-body territories."""
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import gaussian
from skimage.segmentation import expand_labels, watershed

T_CELL, TUMOUR = "T", "Tumour"


def smooth(img, voxel, sigma_um):
    return gaussian(img, sigma=[sigma_um / v for v in voxel], preserve_range=True)


def ball(radius_um, voxel):
    """Ellipsoidal footprint that is a sphere of `radius_um` in physical space."""
    r = [max(1, int(round(radius_um / v))) for v in voxel]
    grids = np.ogrid[tuple(slice(-n, n + 1) for n in r)]
    return sum((g / n) ** 2 for g, n in zip(grids, r)) <= 1


def background(img, cells=None, voxel=None, margin_um=2.0):
    """(median, robust SD) of the signal outside the cells; whole image if no mask is given."""
    out = np.ones(img.shape, bool)
    if cells is not None:
        out = ~ndi.binary_dilation(cells > 0, ball(margin_um, voxel))
    if out.sum() < 100:
        out = np.ones(img.shape, bool)
    values = img[out]
    return float(np.median(values)), float(1.4826 * np.median(np.abs(values - np.median(values))))


def positive(img, voxel, snr=3.0, sigma_um=0.3):
    """Mask of voxels clearly above the (antibody) background: bg + snr * robust SD.

    The background level and its noise are taken from the *raw* image - smoothing shrinks the
    noise, so a threshold derived from the smoothed image would sit far too close to the
    background and let the whole field through.
    """
    bg, sd = background(img)
    return smooth(img, voxel, sigma_um) > bg + snr * sd


def half_max(img, voxel, snr=3.0, sigma_um=0.3):
    """Like `positive`, but the edge is put half-way between background and the stained signal.

    A fluorescent edge is blurred by the point spread function, so the background + snr * SD
    level of `positive` sits a few hundred nm *outside* the real boundary and inflates every
    cell. The half-maximum level is where a blurred step actually crosses its own mid-point.
    """
    sm = smooth(img, voxel, sigma_um)
    bg, _ = background(img)
    inside = positive(img, voxel, snr)
    if not inside.any():
        return inside
    return sm > 0.5 * (bg + np.median(sm[inside]))


def classify(nuclei, gfp, voxel, snr=3.0, gfp_fraction=0.3, shell_um=1.0):
    """{label: ('T'|'Tumour', GFP+ fraction)}; T cells are the DAPI+ GFP+ nuclei."""
    pos = positive(gfp, voxel, snr)
    footprint = ball(shell_um, voxel)
    out = {}
    for label, sl in slices(nuclei):
        region = ndi.binary_dilation(nuclei[sl] == label, footprint)
        frac = float(pos[sl][region].mean()) if region.any() else 0.0
        out[label] = (T_CELL if frac >= gfp_fraction else TUMOUR, frac)
    return out


def cell_bodies(nuclei, gfp, mcherry, dapi, voxel, snr=3.0, open_um=0.4, expand_um=0.0,
                split="shape"):
    """Grow every nucleus into a cell body and split the touching cells apart.

    The footprint of a cell is taken from **every channel at once** - DAPI, GFP and mCherry
    thresholded at their half maximum, closed, hole-filled, and cleaned of isolated speckle.
    Using one channel per cell type fails whenever that channel is weak: a tumour cell with a
    dim antigen stain then collapses onto its nucleus and the T cell next to it takes the space
    in between, which turns a wide apposition into a sliver.

    Two cells that touch share one blob of foreground, and `split` decides where the border goes:

    - "shape" cuts it at the waist of that blob (watershed on the distance transform), which is
      where two cells pressed together actually meet, whatever their relative size. Default.
    - "nucleus" puts it half-way between the two nuclei. Stable, but it misplaces the border
      between cells of unequal size, pushing it into the larger one.

    Returns (bodies, cores): the territories grown by `expand_um` and the ungrown ones. Contacts
    are looked for in `bodies`; intensities are read around `cores`.
    """
    fg = np.zeros(dapi.shape, bool)
    for channel in (dapi, gfp, mcherry):  # filled per channel: a hollow surface stain
        # only becomes a body once its own shell is closed and filled
        fg |= fill(ndi.binary_closing(half_max(channel, voxel, snr), ball(0.6, voxel)))
    if open_um > 0:
        fg = ndi.binary_opening(fg, ball(open_um, voxel))  # isolated noise speckle
    fg |= nuclei > 0
    blobs, _ = ndi.label(fg)
    keep = np.unique(blobs[nuclei > 0])
    fg = np.isin(blobs, keep[keep > 0])

    if split == "shape":
        elevation = -smooth(ndi.distance_transform_edt(fg, sampling=voxel), voxel, 0.3)
    else:
        elevation = ndi.distance_transform_edt(nuclei == 0, sampling=voxel)
    cores = watershed(elevation, nuclei, mask=fg)
    bodies = expand_labels(cores, distance=expand_um, spacing=voxel) if expand_um > 0 else cores
    return bodies, cores


def fill(mask):
    """Fill enclosed holes in 3D and, additionally, within each z plane.

    A surface stain is a hollow shell; filling it turns it into a cell body. With a coarse z step
    the shell is broken between planes and the 3D fill leaks, so each plane is filled on its own
    as well - in a single plane the shell is still a closed ring.
    """
    filled = ndi.binary_fill_holes(mask)
    for z, plane in enumerate(mask):
        filled[z] |= ndi.binary_fill_holes(plane)
    return filled


def clipped(mask, axes="zyx"):
    """True if the cell runs into the edge of the field of view along any of `axes`.

    "xy" asks only about the sides of the field, "zyx" also about the first and last z plane.
    A stack that is barely deeper than a cell clips most cells in z, so excluding those would
    throw away the whole experiment; they are flagged instead.
    """
    for axis, name in enumerate("zyx"):
        if name in axes and (mask.take(0, axis).any() or mask.take(-1, axis).any()):
            return True
    return False


def slices(labels):
    """(label, bounding-box slice) for every label present."""
    for label, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is not None:
            yield label, sl
