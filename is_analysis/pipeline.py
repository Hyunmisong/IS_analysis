"""End-to-end analysis: *_stack z-stacks -> synapse table, QC figures, group statistics and plot."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import stats

from . import plotting, rois, segmentation, synapse
from .io import list_stacks, load_stack

KEY_METRIC = "Contact_Area_um2"
MIN_PLANES = 4


def analyse_image(path, out, voxel_override=None, channels=(0, 1, 2), opts=None):
    """Segment one image, measure every synapse in it and write its figures."""
    name = Path(path).stem
    stack, voxel = load_stack(path)
    if voxel_override:
        voxel = voxel_override
    if stack.shape[1] < MIN_PLANES:  # overview tiles and single planes are not z-stacks
        print("{}: skipped, only {} z plane(s)".format(Path(path).name, stack.shape[1]))
        return None
    if max(channels) >= stack.shape[0]:
        raise SystemExit("{} has {} channel(s); --channels {} is out of range".format(
            Path(path).name, stack.shape[0], ",".join(str(c) for c in channels)))
    stack = np.stack([stack[c] for c in channels])
    dapi, gfp, mcherry = stack

    nuclei = segmentation.segment_nuclei(dapi, voxel, opts["min_nucleus_um3"],
                                         opts["min_separation_um"])
    classes = segmentation.classify(nuclei, gfp, voxel, opts["snr"], opts["gfp_fraction"])
    bodies, cores = segmentation.cell_bodies(nuclei, gfp, mcherry, dapi, voxel, opts["snr"],
                                             opts["open_um"], opts["expand_um"], opts["split"])
    clipped_z = {}
    for label, _ in list(segmentation.slices(bodies)):
        mask = bodies == label
        clipped_z[label] = segmentation.clipped(mask, "z")
        if opts["border"] != "none" and segmentation.clipped(mask, opts["border"]):
            bodies[mask] = 0
            cores[mask] = 0
            nuclei[nuclei == label] = 0
    classes = {k: v for k, v in classes.items() if (bodies == k).any()}

    rows, synapses = synapse.measure(stack, nuclei, bodies, cores, classes, voxel,
                                     opts["min_area_um2"], opts["shell_um"])
    n_t = sum(1 for kind, _ in classes.values() if kind == segmentation.T_CELL)
    for r in rows:
        r.update(Image=name, Voxel_Z_um=voxel[0], Voxel_XY_um=voxel[2],
                 N_T_Cells=n_t, N_Tumour_Cells=len(classes) - n_t,
                 Clipped_Z=clipped_z[r["T_Cell"]] or clipped_z[r["Tumour_Cell"]])

    _write_masks(out / "masks", name, bodies, nuclei, synapses, voxel)
    rois.write(out / "rois" / f"{name}.zip", bodies, synapses, classes)
    plotting.plot_image_qc(
        out / "segmentation_qc" / f"{name}.png", stack, bodies, classes, rows, voxel,
        f"{name}  |  T cells: {n_t}  tumour cells: {len(classes) - n_t}  synapses: {len(rows)}")
    for r in rows:
        plotting.plot_synapse(
            out / "per_synapse" / f"{name}_IS{r['IS']:02d}.png", stack, bodies, r, voxel,
            "{}  IS {}  T{}-Tumour{}  area={:.2f} um2  d={:.2f} um".format(
                name, r["IS"], r["T_Cell"], r["Tumour_Cell"], r[KEY_METRIC],
                r["Contact_Diameter_um"]))
    return rows


def _write_masks(out, name, bodies, nuclei, synapses, voxel):
    """Label images as ImageJ-readable 16-bit stacks, so the result can be checked in Fiji."""
    out.mkdir(parents=True, exist_ok=True)
    for tag, arr in (("bodies", bodies), ("nuclei", nuclei), ("synapses", synapses)):
        tifffile.imwrite(out / f"{name}_{tag}.tif", arr.astype(np.uint16), imagej=True,
                         compression="lzw",  # label images are ~100x smaller compressed
                         resolution=(1 / voxel[2], 1 / voxel[1]),
                         metadata={"spacing": voxel[0], "unit": "um", "axes": "ZYX"})


def group_stats(df, value, control):
    """Per condition: n, median, and Mann-Whitney U against the control condition."""
    groups = [g[value].dropna() for _, g in df.groupby("Condition")]
    kw_p = stats.kruskal(*groups).pvalue if len(groups) > 1 else float("nan")
    ctrl = df.loc[df.Condition == control, value].dropna()
    rows = []
    for condition, g in df.groupby("Condition"):
        v = g[value].dropna()
        p = float("nan")
        if condition != control and len(v) and len(ctrl):
            p = stats.mannwhitneyu(v, ctrl).pvalue
        rows.append(dict(Condition=condition, n_synapses=len(v), n_images=g.Image.nunique(),
                         median=v.median(), mean=v.mean(), p_vs_control=p))
    return pd.DataFrame(rows), kw_p


def _conditions(df, samples_path):
    """Attach Condition from metadata/samples.csv; fall back to the filename prefix."""
    path = Path(samples_path)
    if path.exists():
        samples = pd.read_csv(path)
        merged = df.merge(samples, left_on="Image", right_on="image", how="left")
        merged["Condition"] = merged["condition"].fillna("unmapped")
        return merged.drop(columns=[c for c in ("image", "condition") if c in merged])
    df["Condition"] = df.Image.str.split(r"[-_]", regex=True).str[0]
    return df


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="data/raw", help="folder with the image files")
    ap.add_argument("--pattern", default="*_stack*", help="filename glob, without the extension")
    ap.add_argument("--samples", default="metadata/samples.csv", help="image -> condition table")
    ap.add_argument("--out", default="results")
    ap.add_argument("--channels", default="0,1,2", help="0-based DAPI,GFP,mCherry channel indices")
    ap.add_argument("--voxel", default="", help="override dz,dy,dx in um (default: file metadata)")
    ap.add_argument("--control", default="", help="condition used as reference in the statistics")
    ap.add_argument("--metric", default=KEY_METRIC, help="metric that is plotted and tested")
    ap.add_argument("--snr", type=float, default=3.0,
                    help="GFP/mCherry positivity cut-off in robust SDs above background")
    ap.add_argument("--gfp-fraction", type=float, default=0.3,
                    help="GFP+ voxel fraction around a nucleus needed to call it a T cell")
    ap.add_argument("--min-nucleus-um3", type=float, default=20.0, help="smallest accepted nucleus")
    ap.add_argument("--min-separation-um", type=float, default=4.0,
                    help="smallest distance between two nuclear centres; raise it if one nucleus "
                         "is split in two, lower it if two nuclei are merged")
    ap.add_argument("--open-um", type=float, default=0.4,
                    help="isolated specks smaller than this are removed from the cell footprint")
    ap.add_argument("--split", choices=("shape", "nucleus"), default="shape",
                    help="put the border between two touching cells at the waist of their "
                         "shared footprint (shape) or half-way between their nuclei (nucleus)")
    ap.add_argument("--expand-um", type=float, default=0.0,
                    help="extra territory growth, if real contacts are missed (inflates the area)")
    ap.add_argument("--min-area-um2", type=float, default=0.5,
                    help="contacts smaller than this are not counted as a synapse")
    ap.add_argument("--shell-um", type=float, default=0.5,
                    help="half-thickness of the shell used for the mCherry readout")
    ap.add_argument("--border", choices=("xy", "zyx", "none"), default="xy",
                    help="drop cells clipped by the sides of the field (xy), by any edge "
                         "including the first/last z plane (zyx), or keep everything (none)")
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    files = list_stacks(a.data_dir, a.pattern)
    if not files:
        raise SystemExit("no image matching '{}' in {}".format(a.pattern, a.data_dir))
    channels = tuple(int(c) for c in a.channels.split(","))
    voxel = tuple(float(v) for v in a.voxel.split(",")) if a.voxel else None
    opts = {k: getattr(a, k) for k in
            ("snr", "gfp_fraction", "min_nucleus_um3", "min_separation_um",
             "expand_um", "min_area_um2", "shell_um", "border", "open_um", "split")}

    rows = []
    for path in files:
        found = analyse_image(path, out, voxel, channels, opts)
        if found is not None:
            print("{}: {} synapse(s)".format(path.name, len(found)), flush=True)
        rows += found or []
    if not rows:
        raise SystemExit("no T cell / tumour cell contact found. Check results/segmentation_qc/, "
                         "then try a lower --gfp-fraction or --snr, or --expand-um 0.3")

    df = _conditions(pd.DataFrame(rows), a.samples)
    front = ["Image", "Condition", "IS", "T_Cell", "Tumour_Cell", KEY_METRIC, "Contact_Diameter_um",
             "Contact_Fraction_T", "mCherry_Enrichment"]
    df = df[front + [c for c in df.columns if c not in front]]
    df.to_csv(out / "synapse_table.csv", index=False)
    try:
        df.to_excel(out / "synapse_table.xlsx", index=False)
    except ImportError:
        pass

    control = a.control or df.Condition.iloc[0]
    stats_df, kw_p = group_stats(df, a.metric, control)
    stats_df.to_csv(out / "group_stats.csv", index=False)
    plotting.plot_groups(
        out / "synapse_size.png", df, a.metric, control,
        dict(zip(stats_df.Condition, stats_df.p_vs_control)),
        "n = {} synapses, Kruskal-Wallis p = {:.3g}\nvs control '{}' "
        "(Mann-Whitney, see the caveat in the README)".format(len(df), kw_p, control))
    print(stats_df.to_string(index=False))
    print("\nwrote {} and {}".format(out / "synapse_table.csv", out / "synapse_size.png"))
