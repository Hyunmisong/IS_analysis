"""Loading of multi-channel z-stacks (.czi from the microscope, .tif exported by ImageJ/Fiji)."""
import re
from pathlib import Path

import numpy as np
import tifffile


def list_stacks(data_dir, pattern="*_stack*"):
    """Image files matching `pattern` (default: the `*_stack*` exports), .czi and .tif."""
    d = Path(data_dir)
    files = [p for ext in ("czi", "tif", "tiff") for p in d.glob(f"{pattern}.{ext}")]
    return sorted(set(files), key=lambda p: p.name)


def load_stack(path, n_channels=3):
    """Return (stack[C, Z, Y, X] float32, voxel size (dz, dy, dx) in um)."""
    path = Path(path)
    if path.suffix.lower() == ".czi":
        return _load_czi(path)
    return _load_tif(path, n_channels)


def _order(arr, axes, n_channels):
    """Drop singleton/unused axes and transpose to C, Z, Y, X."""
    arr = arr[tuple(0 if a not in "CZYX" else slice(None) for a in axes)]
    axes = "".join(a for a in axes if a in "CZYX")
    for a in "CZ":  # restore axes the file does not carry
        if a not in axes:
            arr, axes = arr[None], a + axes
    arr = np.transpose(arr, [axes.index(a) for a in "CZYX"])
    if arr.shape[0] == 1 and arr.shape[1] % n_channels == 0 and arr.shape[1] > n_channels:
        # ImageJ wrote the hyperstack as one flat plane list (channel changes fastest)
        arr = arr[0].reshape(-1, n_channels, *arr.shape[2:]).transpose(1, 0, 2, 3)
    return np.ascontiguousarray(arr, dtype=np.float32)


def _load_czi(path):
    import czifile

    czi = czifile.CziFile(str(path))
    stack = _order(czi.asarray(), czi.axes, 3)
    meta = czi.metadata()
    um = {}
    for axis, value in re.findall(r'<Distance Id="([XYZ])">\s*<Value>([^<]+)</Value>', meta):
        um[axis] = float(value) * 1e6
    voxel = (um.get("Z", um.get("X", 1.0)), um.get("Y", 1.0), um.get("X", 1.0))
    return stack, voxel


def _load_tif(path, n_channels):
    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        stack = _order(series.asarray(), series.axes, n_channels)
        ij = tif.imagej_metadata or {}
        dz = float(ij.get("spacing", 1.0))
        tags = tif.pages[0].tags
        dxy = 1.0
        if "XResolution" in tags:
            num, den = tags["XResolution"].value
            dxy = den / num if num else 1.0
    return stack, (dz, dxy, dxy)

