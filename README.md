# Immunological synapse analysis (DAPI / GFP / mCherry confocal z-stacks)

A pipeline that finds every **T cell – tumour cell contact** in a confocal z-stack and measures
**how large that contact is**, together with how much of the **mCherry-stained antigen (BCMA)**
sits in it.

| cell type | how it is recognised |
|---|---|
| **Tumour cell** | DAPI+, GFP− ; its body is outlined by the mCherry antigen stain on its surface |
| **T cell** | DAPI+ **and** GFP+ ; GFP fills the cytosol, so the GFP+ volume *is* the T cell |
| **Immunological synapse (IS)** | the surface where a T cell body and a tumour cell body touch |

![example synapse](docs/example_synapse.png)

## Why Python and not an ImageJ macro

The measurement is a **3D** one: splitting touching nuclei, growing each nucleus into a cell body,
and computing the area of the surface shared by two of those bodies. In the ImageJ macro language
that means leaning on 3D plugins (3D ImageJ Suite, MorphoLibJ) and still writing the contact-area
step yourself, with no straightforward way to script the per-condition statistics afterwards.
In Python, `scikit-image` and `scipy.ndimage` do the 3D work directly, and the whole study reruns
with one command when a parameter changes.

**ImageJ still does the two things it is best at**, and the repository ships macros for both:

| Step | Tool |
|---|---|
| Convert the microscope files into the `*_stack.tif` hyperstacks | [`tools/export_stacks.ijm`](tools/export_stacks.ijm) (Fiji, Bio-Formats) |
| Measure | `python scripts/run_analysis.py` |
| Look at the segmentation plane by plane | [`tools/check_masks.ijm`](tools/check_masks.ijm) (Fiji) |

The pipeline also reads `.czi` directly, so the Fiji export step is optional.
Every mask it produces is written as an ImageJ-readable 16-bit stack in `results/masks/`.

## Usage

```bash
pip install -r requirements.txt
# put the *_stack.tif (or .czi) files in data/raw/, then
python scripts/run_analysis.py --data-dir data/raw
```

The image → condition mapping lives in [`metadata/samples.csv`](metadata/samples.csv); edit it for
your own experiment. Without it (`--samples none`) the condition is taken from the part of the
filename before the first `-` or `_`.

Useful options:

| Option | Default | Meaning |
|---|---|---|
| `--pattern` | `*_stack*` | which files to analyse (without the extension) |
| `--channels` | `0,1,2` | 0-based indices of DAPI, GFP, mCherry **in the order they were acquired** |
| `--voxel` | from file metadata | override `dz,dy,dx` in µm |
| `--snr` | `3.0` | how many robust SDs above background a voxel must be to count as stained |
| `--gfp-fraction` | `0.3` | GFP+ fraction around a nucleus needed to call it a T cell |
| `--min-separation-um` | `3.0` | smallest distance between two nuclear centres (splitting touching nuclei) |
| `--expand-um` | `0.0` | extra growth of each cell territory; raise to `0.3` if real contacts are missed |
| `--min-area-um2` | `0.5` | contacts smaller than this are not counted as a synapse |
| `--metric` | `Contact_Area_um2` | which column is plotted and tested |
| `--control` | first condition | reference group of the statistics |

## Pipeline

| Step | What it does | Code |
|---|---|---|
| 1 | Load the stack (`.czi` or ImageJ `.tif`), read the voxel size from the metadata | `io.py` |
| 2 | **Nuclei (DAPI)**: Gaussian smoothing → Otsu → hole filling → touching nuclei split by a watershed on the 3D distance map (anisotropic voxels handled throughout) | `segmentation.py` |
| 3 | **Classify**: fraction of GFP+ voxels in a 1 µm shell around each nucleus → T cell or tumour cell | `segmentation.py` |
| 4 | **Cell bodies**: the T cell body is the GFP+ volume, the tumour body is the filled mCherry+ surface; each type is split among its own nuclei. Voxels claimed by both channels go to the nearer nucleus | `segmentation.py` |
| 5 | **Contacts**: every pair of touching bodies; T–T and tumour–tumour pairs are discarded | `synapse.py` |
| 6 | **Size**: the contact area from a marching-cubes mesh of the T cell, taking the triangles that face the tumour cell | `synapse.py` |
| 7 | **Antigen**: mCherry in a ±0.5 µm band on the tumour surface, at the synapse vs. over the rest of that surface, both background-subtracted | `synapse.py` |
| 8 | Per-image and per-synapse figures, pooled table, group statistics, final plot | `plotting.py`, `pipeline.py` |

Cells clipped by the edge of the field of view (x, y or the first/last z plane) are dropped, since
their volume and surface are not measurable; `--keep-border` keeps them.

### Where the boundary of a cell is put

A fluorescent edge is blurred by the microscope, so a "background + 3 SD" threshold sits a few
hundred nm *outside* the real boundary and inflates every cell. The body masks therefore use the
**half-maximum level** — half-way between the background and the stained signal — which is where a
blurred step actually crosses its own mid-point. On the synthetic test below this brings the cell
volumes to within 10 % of the truth instead of ~80 % too large.

### Background

mCherry is detected with an antibody, so there is a diffuse background. It is estimated **per
image** as the median of the voxels more than 2 µm away from any detected cell, and subtracted
from every mCherry number. `mCherry_SNR` (tumour surface signal / background noise) is reported
per synapse — if it is close to 1, the staining failed in that image and its antigen numbers mean
nothing. The *contact area* does not depend on the mCherry level, only on the segmented boundary.

## Outputs (`results/`)

| File | Description |
|---|---|
| `synapse_table.csv` / `.xlsx` | **one row per T cell – tumour cell contact**, pooled over all images |
| `per_synapse/{image}_IS01.png` | the contact plane of one synapse: merge, mCherry, the two cell bodies, contact outlined in white |
| `segmentation_qc/{image}.png` | projection with every cell outlined (cyan = T cell, yellow = tumour) and every synapse marked |
| `masks/{image}_bodies.tif`, `_nuclei.tif` | the label images, ImageJ-readable, for checking in Fiji |
| `group_stats.csv` | per condition: n, median, mean, Mann-Whitney p vs the control |
| `synapse_size.png` | the final box plot |

![segmentation QC](docs/example_qc.png)

### Columns of `synapse_table.csv`

| Column | Meaning |
|---|---|
| `Contact_Area_um2` | **the IS size**: area of the shared surface, from the marching-cubes mesh |
| `Contact_Diameter_um` | diameter of a circle of that area (`2·√(A/π)`) |
| `Contact_Length_um`, `Contact_Width_um` | the two largest principal extents of the contact patch |
| `Contact_Fraction_T` | contact area ÷ total T cell surface — the size-independent version of the metric |
| `Contact_Area_Voxelface_um2` | the same area counted as voxel faces, as an independent cross-check |
| `T_Volume_um3`, `Tumour_Volume_um3`, `T_Surface_um2` | cell geometry |
| `mCherry_IS`, `mCherry_TumourSurface` | background-subtracted antigen in the surface band, at the synapse and elsewhere |
| `mCherry_Enrichment` | their ratio; > 1 = antigen concentrated at the synapse |
| `mCherry_Background`, `mCherry_SNR` | staining quality of that image |
| `N_T_Cells`, `N_Tumour_Cells` | cells kept in that image, i.e. how crowded the field was |

**`Contact_Fraction_T` is usually the fairer comparison**: a bigger T cell makes a bigger contact
without being any more activated, and this column divides that out.

## Validation

`scripts/make_test_data.py` builds synthetic pairs — two spheres flattened against a common plane,
with a thin gap, a surface antigen rim, a diffuse antibody background and Poisson noise — so the
true contact area is exactly π·r². `scripts/validate.py` runs the whole pipeline on them:

```bash
python scripts/validate.py --dz 0.4
```

| true area (µm²) | measured, dz = 0.4 µm | measured, dz = 1.0 µm | true contact diameter (µm) | measured `Contact_Length_um`, dz = 0.4 / 1.0 |
|---|---|---|---|---|
| 7.07 | 9.90 (×1.40) | 7.66 (×1.08) | 3.0 | 3.2 / 2.8 |
| 19.63 | 22.72 (×1.16) | 22.46 (×1.14) | 5.0 | 4.9 / 4.8 |
| 38.48 | 48.55 (×1.26) | 49.15 (×1.28) | 7.0 | 6.8 / 6.4 |

Cell volumes come out within ~10 % (T cell 268 µm³ true, 256–267 measured; tumour 1437 µm³ true,
1470–1544 measured).

So: **`Contact_Area_um2` is systematically 10–40 % too large in absolute terms** — the segmented
boundary cannot be more precise than the point spread function, and a contact is the small
difference between two large surfaces. The ranking between conditions is preserved, which is what
the comparison needs. `Contact_Length_um` tracks the true contact diameter to within ~0.5 µm and
is the least biased single number here.

Run `python scripts/validate.py` after changing anything in `is_analysis/`.

## Statistical caveat: the p-values are overestimated

The significance marks (Kruskal-Wallis across conditions, Mann-Whitney U vs. the control) treat
**each synapse as an independent observation**. It is not: synapses from the same image, field and
dish share culture, staining and imaging conditions, and several synapses can even involve the
same T cell. This is **pseudoreplication**, so the effective sample size is much smaller than the
number of rows and the **p-values are too small (anti-conservative)**.

Read `*`, `**`, `***` as exploratory. A rigorous test needs independent biological replicates
(dish/well), with either the replicate mean as the unit of analysis or a mixed-effects model with
image/dish as a random effect.

## Known limitations

- **Only contacts that are in the stack are found.** A cell clipped by the field of view is
  dropped, so synapses at the edge are missed and n is smaller than what the eye counts.
- **z resolution.** A 1 µm z step still worked on the synthetic data, but it under-samples the
  ~0.8 µm antigen rim; dz ≤ 0.4 µm is strongly preferred at 63×, and the same dz must be used for
  all conditions that are compared.
- **A contact is not proof of a synapse.** Two cells can touch by chance, and a real IS has
  molecular polarisation this pipeline does not test. The small contacts in particular (close to
  `--min-area-um2`) are better read as "cells in contact" than as mature synapses.
- **Tumour bodies depend on the mCherry stain.** Where the antibody signal is weak the tumour body
  falls back to a 1.5 µm shell around the nucleus and its contact is then underestimated. Check
  `mCherry_SNR` and the figures in `segmentation_qc/` before trusting an image.
- **GFP classification is a threshold.** A dim GFP+ T cell is called a tumour cell; check
  `segmentation_qc/` (cyan vs. yellow outlines) and tune `--gfp-fraction` if the colours are wrong.
- The raw microscopy files are not tracked in git (~40 MB each, see `.gitignore`).
