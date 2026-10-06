"""Write synthetic *_stack.tif images with a known contact area, to validate the pipeline.

Each image holds one T cell / tumour cell pair. Both cells are spheres flattened against a common
plane, as two cells pressed together are, so the contact is a flat disc of a radius we choose and
the true contact area is exactly pi * r^2. A thin gap is left between the two flat faces, the
antigen (mCherry) sits on the tumour surface on top of a diffuse antibody background, and every
channel gets Poisson noise at a realistic photon count.

    python scripts/make_test_data.py --out data/test --dz 0.4
    python scripts/run_analysis.py --data-dir data/test --samples none --out results_test
"""
import argparse
from pathlib import Path

import numpy as np
import tifffile

VOXEL_XY = 0.2
RNG = np.random.default_rng(0)


def _grid(shape):
    return np.ogrid[tuple(slice(0, n) for n in shape)]


def sphere(shape, centre, radius, voxel):
    return sum(((g - c) * v / radius) ** 2
               for g, c, v in zip(_grid(shape), centre, voxel)) <= 1.0


def flattened(shape, centre, radius, voxel, plane_x, side):
    """Sphere cut off by the plane x = plane_x; `side` = -1 keeps the lower-x half."""
    x = _grid(shape)[2] * voxel[2]
    return sphere(shape, centre, radius, voxel) & ((x - plane_x) * side >= 0)


def make_pair(contact_r=2.5, r_t=4.0, r_tumour=7.0, gap_um=0.3, voxel=(0.4, VOXEL_XY, VOXEL_XY)):
    """(stack[C, Z, Y, X] uint16, true contact area in um^2)."""
    h_t, h_tumour = np.sqrt(r_t**2 - contact_r**2), np.sqrt(r_tumour**2 - contact_r**2)
    size_um = np.array([2 * r_tumour + 6, 2 * r_tumour + 6, 2 * (r_t + r_tumour) + 8])
    shape = tuple(int(round(s / v)) for s, v in zip(size_um, voxel))
    cy, cz = np.array(shape[:2]) / 2
    plane = size_um[2] / 2
    c_t = (cy, cz, (plane - gap_um / 2 - h_t) / voxel[2])
    c_tumour = (cy, cz, (plane + gap_um / 2 + h_tumour) / voxel[2])

    t_body = flattened(shape, c_t, r_t, voxel, plane - gap_um / 2, -1)
    tumour_body = flattened(shape, c_tumour, r_tumour, voxel, plane + gap_um / 2, +1)
    dapi = (sphere(shape, c_t, r_t * 0.5, voxel) * 2200.0
            + sphere(shape, c_tumour, r_tumour * 0.5, voxel) * 2600.0)
    gfp = t_body * 1600.0
    rim = tumour_body & ~_shrink(tumour_body, shape, c_tumour, r_tumour, voxel,
                                 plane + gap_um / 2, 0.8)
    mcherry = rim * 1400.0 + 250.0  # antigen on the tumour surface + antibody background

    stack = np.stack([dapi, gfp, mcherry]) + 120.0
    stack = RNG.poisson(stack / 6.0) * 6.0  # shot noise at a realistic photon count
    return stack.astype(np.uint16), float(np.pi * contact_r**2)


def _shrink(body, shape, centre, radius, voxel, plane_x, by_um):
    """The same flattened sphere, `by_um` smaller in every direction (gives a surface rim)."""
    return flattened(shape, centre, radius - by_um, voxel, plane_x + by_um, +1)


def write_images(out, dz, radii):
    """Write one *_stack.tif per contact radius; returns the true contact areas in um^2."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    areas = []
    for i, r in enumerate(radii, start=1):
        stack, area = make_pair(contact_r=r, voxel=(dz, VOXEL_XY, VOXEL_XY))
        name = out / f"{i}-1_stack.tif"
        tifffile.imwrite(name, stack.transpose(1, 0, 2, 3), imagej=True,
                         resolution=(1 / VOXEL_XY, 1 / VOXEL_XY),
                         metadata={"spacing": dz, "unit": "um", "axes": "ZCYX"})
        print(f"{name.name}: contact radius {r} um -> true contact area = {area:.2f} um^2")
        areas.append(area)
    return areas


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/test")
    ap.add_argument("--dz", type=float, default=0.4, help="z step in um")
    ap.add_argument("--radii", default="1.5,2.5,3.5", help="contact radii in um, one per image")
    a = ap.parse_args()
    print(f"voxel dz,dy,dx = {(a.dz, VOXEL_XY, VOXEL_XY)} um")
    write_images(a.out, a.dz, [float(v) for v in a.radii.split(",")])


if __name__ == "__main__":
    main()
