"""Detection and measurement of the T cell / tumour cell contact (immunological synapse)."""
import numpy as np
from scipy import ndimage as ndi
from skimage.measure import marching_cubes

from .segmentation import T_CELL, TUMOUR, ball, background


def contact_pairs(bodies, voxel):
    """{(label_a, label_b): staircase contact area in um^2} for every pair of touching cells.

    One pass over the three face directions: two neighbouring voxels that belong to different
    cells share a face, and the face area depends on the direction (dy*dx for a z-face, ...).
    """
    area = {}
    for axis in range(3):
        face = np.prod(voxel) / voxel[axis]
        lo = np.take(bodies, np.arange(bodies.shape[axis] - 1), axis=axis)
        hi = np.take(bodies, np.arange(1, bodies.shape[axis]), axis=axis)
        touch = (lo > 0) & (hi > 0) & (lo != hi)
        if not touch.any():
            continue
        a, b = lo[touch], hi[touch]
        pairs, counts = np.unique(np.stack([np.minimum(a, b), np.maximum(a, b)]), axis=1,
                                  return_counts=True)
        for (p, q), n in zip(pairs.T, counts):
            key = (int(p), int(q))
            area[key] = area.get(key, 0.0) + n * face
    return area


def interface_mask(bodies, a, b, voxel):
    """Voxels of cell `a` that touch cell `b`, plus the mirror voxels of `b` (the contact shell)."""
    ma, mb = bodies == a, bodies == b
    cross = ndi.generate_binary_structure(3, 1)  # one voxel step along each axis
    return (ma & ndi.binary_dilation(mb, cross)) | (mb & ndi.binary_dilation(ma, cross))


def surface_area(mask, voxel):
    """Area of the mask surface in um^2 from a marching-cubes mesh (no staircase bias)."""
    padded = np.pad(mask.astype(np.float32), 1)
    try:
        verts, faces, *_ = marching_cubes(padded, level=0.5, spacing=voxel)
    except (RuntimeError, ValueError):
        return float("nan")
    tri = verts[faces]
    return float(0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]),
                                      axis=1).sum())


def contact_area_mesh(bodies, a, b, voxel):
    """Contact area from the mesh of cell `a`: triangles that face cell `b`.

    Less biased than counting voxel faces, which overestimates an oblique surface by up to ~1.5x.
    """
    ma = bodies == a
    padded = np.pad(ma.astype(np.float32), 1)
    try:
        verts, faces, *_ = marching_cubes(padded, level=0.5, spacing=voxel)
    except (RuntimeError, ValueError):
        return float("nan"), float("nan")
    tri = verts[faces]
    areas = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    centre = tri.mean(axis=1) / np.array(voxel) - 1.0  # mesh coords -> voxel index of `bodies`
    near_b = ndi.binary_dilation(bodies == b, ndi.generate_binary_structure(3, 1))
    idx = np.clip(np.round(centre).astype(int).T, 0, np.array(bodies.shape)[:, None] - 1)
    at_contact = near_b[tuple(idx)]
    return float(areas[at_contact].sum()), float(areas.sum())


def extents_um(mask, voxel):
    """The two largest principal extents of the contact patch in um (its 'width' and 'length')."""
    pts = np.argwhere(mask) * np.array(voxel)
    if len(pts) < 3:
        return float("nan"), float("nan")
    pts = pts - pts.mean(0)
    proj = pts @ np.linalg.svd(pts, full_matrices=False)[2].T
    spread = sorted(np.ptp(proj, axis=0), reverse=True)
    return float(spread[0]), float(spread[1])


def antigen_at_synapse(mcherry, bodies, tumour, contact, voxel, bg, bg_sd, shell_um=0.5):
    """mCherry (BCMA) in the tumour surface band at the synapse vs over the rest of that band.

    The band straddles the tumour boundary by +/- `shell_um`, so it captures the stain whether
    the antibody signal sits just inside or just outside the segmented edge. Both numbers are
    background-subtracted, so `Enrichment` > 1 means the antigen is concentrated at the synapse.
    """
    footprint = ball(shell_um, voxel)
    tumour_mask = bodies == tumour
    band = ndi.binary_dilation(tumour_mask, footprint) & ~ndi.binary_erosion(tumour_mask, footprint)
    near_is = ndi.binary_dilation(contact, footprint)
    at_is, elsewhere = band & near_is, band & ~near_is
    is_mean = float(mcherry[at_is].mean()) - bg if at_is.any() else float("nan")
    rest = float(mcherry[elsewhere].mean()) - bg if elsewhere.any() else float("nan")
    whole = float(mcherry[band].mean()) - bg if band.any() else float("nan")
    return dict(mCherry_IS=is_mean, mCherry_TumourSurface=rest, mCherry_TumourMean=whole,
                mCherry_Enrichment=is_mean / rest if rest and rest > 0 else float("nan"),
                mCherry_Background=bg, mCherry_SNR=whole / bg_sd if bg_sd else float("nan"))


def _crop(bodies, a, b, pad=2):
    """Bounding box that contains both cells completely, so areas and volumes stay exact."""
    where = np.argwhere((bodies == a) | (bodies == b))
    lo = np.maximum(where.min(0) - pad, 0)
    hi = np.minimum(where.max(0) + pad + 1, bodies.shape)
    return tuple(slice(int(l), int(h)) for l, h in zip(lo, hi))


def measure(stack, nuclei, bodies, cores, classes, voxel, min_area_um2=0.5, shell_um=0.5):
    """(rows, synapse label image): one row per T cell / tumour cell contact, largest first.

    The label image carries the contact voxels of synapse *i* under the value *i*, i.e. the same
    numbering as the `IS` column and the per-synapse figures, so it can be overlaid in Fiji.
    """
    dapi, gfp, mcherry = stack
    bg, bg_sd = background(mcherry, cores, voxel)
    volume = float(np.prod(voxel))
    found = []
    for (a, b), staircase in contact_pairs(bodies, voxel).items():
        kinds = {classes[a][0]: a, classes[b][0]: b}
        if len(kinds) != 2:  # T-T or tumour-tumour contact: not a synapse
            continue
        t, tumour = kinds[T_CELL], kinds[TUMOUR]
        sl = _crop(bodies, t, tumour)
        sub, sub_core, sub_nuc = bodies[sl], cores[sl], nuclei[sl]
        sub_dapi, sub_gfp, sub_mch = dapi[sl], gfp[sl], mcherry[sl]
        contact = interface_mask(sub, t, tumour, voxel)
        area, t_surface = contact_area_mesh(sub, t, tumour, voxel)
        if not (area >= min_area_um2):
            continue
        length, width = extents_um(contact, voxel)
        found.append((sl, contact, dict(
            T_Cell=t, Tumour_Cell=tumour,
            Contact_Area_um2=area,
            Contact_Area_Voxelface_um2=staircase,
            Contact_Diameter_um=2 * np.sqrt(area / np.pi),
            Contact_Length_um=length, Contact_Width_um=width,
            Contact_Fraction_T=area / t_surface if t_surface else float("nan"),
            T_Surface_um2=t_surface,
            T_Volume_um3=float((sub == t).sum()) * volume,
            Tumour_Volume_um3=float((sub == tumour).sum()) * volume,
            GFP_IS=float(sub_gfp[contact].mean()),
            DAPI_T=float(sub_dapi[sub_nuc == t].mean()) if (sub_nuc == t).any() else float("nan"),
            Z_Contact=float(np.argwhere(contact)[:, 0].mean()) + sl[0].start,
            **antigen_at_synapse(sub_mch, sub_core, tumour, contact, voxel, bg, bg_sd, shell_um),
        )))
    found.sort(key=lambda f: -f[2]["Contact_Area_um2"])
    synapses = np.zeros(bodies.shape, np.uint16)
    rows = []
    for index, (sl, contact, row) in enumerate(found, start=1):
        synapses[sl][contact] = index
        rows.append(dict(IS=index, **row))
    return rows, synapses
