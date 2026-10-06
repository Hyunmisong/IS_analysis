"""Validate the pipeline on synthetic pairs whose contact area is known exactly.

Builds the images of `make_test_data.py` in a temporary folder, runs the full analysis on them
and prints measured vs. true contact area. Run it after changing anything in `is_analysis/`.

    python scripts/validate.py              # z step 0.4 um
    python scripts/validate.py --dz 1.0     # the z step of a quick acquisition
"""
import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from is_analysis.pipeline import main as run_analysis  # noqa: E402

import make_test_data  # noqa: E402  (same folder)

RADII = [1.5, 2.5, 3.5]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dz", type=float, default=0.4, help="z step in um")
    ap.add_argument("--keep", default="", help="also copy the figures into this folder")
    a = ap.parse_args()

    work = Path(tempfile.mkdtemp(prefix="is_validate_"))
    raw, out = work / "raw", work / "out"
    true_areas = make_test_data.write_images(raw, a.dz, RADII)
    run_analysis(["--data-dir", str(raw), "--samples", "none", "--out", str(out)])

    df = pd.read_csv(out / "synapse_table.csv").sort_values("Image")
    df["True_Area_um2"] = true_areas
    df["True_Diameter_um"] = [2 * r for r in RADII]
    df["Area_ratio"] = df.Contact_Area_um2 / df.True_Area_um2
    columns = ["Image", "True_Area_um2", "Contact_Area_um2", "Area_ratio",
               "True_Diameter_um", "Contact_Length_um", "T_Volume_um3", "Tumour_Volume_um3"]
    print("\nz step {} um - true T cell volume 268 um3, tumour 1437 um3\n".format(a.dz))
    print(df[columns].round(2).to_string(index=False))

    if a.keep:
        shutil.copytree(out, Path(a.keep), dirs_exist_ok=True)
        print("\nfigures copied to", a.keep)
    else:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
