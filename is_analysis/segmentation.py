"""3D nuclei segmentation, T cell / tumour cell classification and cell-body territories."""
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import remove_small_objects
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


def segment_nuclei(dapi, voxel, min_volume_um3=20.0, min_sep_um=3.0, sigma_um=0.4):
    """Label nuclei in 3D; touching nuclei are split by a distance-transform watershed."""
    sm = smooth(dapi, voxel, sigma_um)
    mask = sm > threshold_otsu(sm)
    mask = ndi.binary_closing(mask, ball(0.5, voxel))
    mask = ndi.binary_fill_holes(mask)
    mask = remove_small_objects(mask, int(min_volume_um3 / np.prod(voxel)))
    dist = ndi.distance_transform_edt(mask, sampling=voxel)
    seeds = ndi.label(_peaks(smooth(dist, voxel, 0.5), mask, voxel, min_sep_um))[0]
    labels = watershed(-dist, seeds, mask=mask)
    return _drop_small(labels, voxel, min_volume_um3)


def _peaks(dist, mask, voxel, min_sep_um):
    """One seed per local maximum of the distance map, maxima closer than min_sep_um merged."""
    footprint = ball(min_sep_um / 2, voxel)
    return mask & (dist >= ndi.maximum_filter(dist, footprint=footprint)) & (dist > 0)


def _drop_small(labels, voxel, min_volume_um3):
    keep = np.bincount(labels.ravel()) * np.prod(voxel) >= min_volume_um3
    keep[0] = False
    out = np.zeros_like(labels)
    for new, old in enumerate(np.flatnonzero(keep), start=1):
        out[labels == old] = new
    return out


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


def cell_bodies(nuclei, classes, gfp, mcherry, voxel, snr=3.0, min_radius_um=1.5,
                expand_um=0.0):
    """Grow every nucleus into a cell body, using the channel that actually marks that cell.

    T cells carry GFP in the cytosol, so their body is the GFP+ volume; the tumour cells are
    outlined by the mCherry antigen stain on their surface, so their body is the filled mCherry+
    volume. Each type is then split among its own nuclei (nearest-nucleus watershed) and the
    Returns (bodies, cores): the territories grown by `expand_um`, which closes the dark gap
    between two cells that are in contact, and the ungrown ones. Contacts are looked for in
    `bodies`; intensities are read around `cores`, whose edge still follows the real staining.
    """
    t_nuclei = np.isin(nuclei, [l for l, (kind, _) in classes.items() if kind == T_CELL]) * nuclei
    tumour_nuclei = np.where(t_nuclei > 0, 0, nuclei)
    t_fg = _body_mask(gfp, voxel, snr, t_nuclei, min_radius_um)
    tumour_fg = _body_mask(mcherry, voxel, snr, tumour_nuclei, min_radius_um)

    to_t = ndi.distance_transform_edt(t_nuclei == 0, sampling=voxel)
    to_tumour = ndi.distance_transform_edt(tumour_nuclei == 0, sampling=voxel)
    both = t_fg & tumour_fg  # claimed by both channels -> goes to the nearer nucleus
    t_fg &= ~(both & (to_tumour < to_t))
    tumour_fg &= ~(both & (to_t <= to_tumour))

    elevation = ndi.distance_transform_edt(nuclei == 0, sampling=voxel)
    cores = watershed(elevation, t_nuclei, mask=t_fg)
    cores = np.where(cores > 0, cores, watershed(elevation, tumour_nuclei, mask=tumour_fg))
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


def _body_mask(img, voxel, snr, own_nuclei, min_radius_um):
    """Signal above background, holes filled, restricted to the blobs that contain a nucleus."""
    mask = fill(ndi.binary_closing(half_max(img, voxel, snr), ball(0.5, voxel)))
    blobs, _ = ndi.label(mask)
    keep = np.unique(blobs[own_nuclei > 0])
    mask = np.isin(blobs, keep[keep > 0])
    return mask | ndi.binary_dilation(own_nuclei > 0, ball(min_radius_um, voxel))


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
