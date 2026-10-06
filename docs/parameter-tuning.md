# How the parameters were chosen

A record of what was tried on the 2026-06-04 co-culture set, what was measured, and what was
kept. **The failed attempts are the point of this file**: several reasonable-sounding changes made
the result worse, and without the numbers somebody would try them again.

Two kinds of measurement are used throughout:

- **Synthetic ground truth** — `scripts/validate.py` builds pairs of cells flattened against a
  common plane, so the true contact area is exactly π·r². Reported as measured ÷ true.
- **Real images** — the number of synapses found, the largest contact, and
  `Contact_Fraction_T` (contact ÷ T cell surface). A contact above ~0.25 of the T cell's whole
  surface is not a synapse but two cells wrapped around each other, usually a segmentation error.
  A *rising* largest contact is the signature of two cells merged into one body.

---

## 1. Which cells to drop at the edge of the field

| `--border` | kept cells (one field) | synapses |
|---|---|---|
| `zyx` (drop anything touching an edge) | 2 of 25 | 0 |
| **`xy` (chosen)** | 25 of 25 | 6 |

The stacks are 13–21 planes at dz = 1 µm, barely deeper than one cell, so almost every cell
touches the first or last plane. Dropping them leaves nothing. They are kept and flagged
`Clipped_Z` instead; their volumes and surfaces are truncated, which is why
`Contact_Fraction_T` should be read with that column in view.

## 2. Where to put the edge of a cell

A blurred fluorescent edge crosses its own half-maximum at the true boundary, while
"background + 3 SD" sits several hundred nm outside it.

| threshold | cell volume vs true |
|---|---|
| background + 3·SD | ~80 % too large |
| **half maximum (chosen)** | within ~20 % |

## 3. How to build the cell body

| method | what happens |
|---|---|
| one channel per cell type — GFP for the T cell, mCherry for the tumour | works only while both stains are strong. Here mCherry is near background, so a tumour collapsed onto a 1.5 µm shell around its nucleus and the GFP blob of the T cell beside it took the space between them. One apposition was reported as a 24 µm² sliver and the tumour came out half its real width |
| **all three channels, each hole-filled on its own (chosen)** | the tumour keeps its real width from DAPI even when its own stain fails. The same pair came out at 43 µm² |

Filling **per channel** matters: a hollow surface stain only becomes a solid body if its own
shell is closed first. Filling the union instead leaves the tumour interior empty.

## 4. Where to cut two cells apart

| `--split` | synthetic, measured ÷ true (3 contacts) | verdict |
|---|---|---|
| `valley` (dip in total fluorescence) | ×47, ×17, ×8.6 | broken — a hollow surface stain has a dark interior, so the basin floods the wrong way and the T cell swallows the tumour |
| `nucleus` (half-way between nuclei) | ×4.9, ×2.7, ×2.1 | the border lands in the middle of the pair, not at the membranes; it is pushed into whichever cell is larger |
| **`shape` (waist of the shared outline, chosen)** | ×1.40, ×1.12, ×1.03 | correct whatever the size difference |

## 5. How far to grow each territory (`--expand-um`)

| `--expand-um` | synthetic area ÷ true | `mCherry_SNR` |
|---|---|---|
| **0.0 (chosen)** | ×1.40, ×1.16, ×1.26 | 7.4 |
| 0.3 | ×2.46, ×1.57, ×1.51 | 1.3 |
| 0.5 | ×4.04, ×2.19, ×1.85 | 0.2 |

Growing the territories closes the gap between two cells that are in contact, but it also pushes
the measured surface outwards: the antigen band then sits outside the real stain and the antigen
readout collapses. The half-maximum boundary already closes gaps of ~0.3 µm, so the default is 0.
Raise it only if real contacts are being missed, and expect the areas to inflate.

## 6. Nucleus detection

This is what everything else rests on: one seed per cell is what keeps two cells apart. A missing
seed makes a cell disappear into its neighbour; a duplicated one cuts a cell in half.

### 6a. Threshold or Cellpose

| detector | result |
|---|---|
| threshold + distance watershed | the DAPI mask breaks into pieces at this photon count. 24 "nuclei" in a plane that holds about 15, some real nuclei split across 2–3 labels, others missing entirely |
| **Cellpose `nuclei` model (chosen)** | round, separated objects that match the DAPI blobs |

Over the whole set:

| detector | synapses | largest contact | `Contact_Fraction_T` > 0.25 |
|---|---|---|---|
| threshold | 52 | 305 µm² | 6 % |
| **Cellpose** | 60 | 194 µm² | 0 % |

### 6b. Cellpose needs the channel smoothed first

| smoothing before Cellpose | nuclei found (field of ~25) |
|---|---|
| none | 5 |
| 0.4 µm | 33 |
| 0.8 µm | 28 |
| **1.0 µm (chosen)** | 27 |

The model was trained on images with far better signal than this; raw, it finds almost nothing.
The smoothing is compensating for the acquisition — see the DAPI note in the README.

### 6c. 2D with stitching, not Cellpose's 3D mode

| mode | nuclei | median z-span | time per image |
|---|---|---|---|
| **stitch 0.25 (chosen)** | 27–35 | 4–5 planes | 27 s |
| `do_3D` with anisotropy | 35 | 9 planes | 141 s |

`do_3D` returns ring-shaped and hollow objects on stacks this anisotropic (1 µm in z against
0.2 µm in x and y), and is five times slower.

---

## What was tried and rejected

### Dropping undersized nuclei as debris

Nucleus volumes span a factor of 25, and a fifth of the detections fall below 20 % of the
expected volume. Dropping them as dead cells or debris looks obviously right. It is not.

| cut-off | nuclei | volume spread (p95/p5) | synapses | largest contact |
|---|---|---|---|---|
| none | 384 | 57 | 65 | 194 µm² |
| **20 µm³ — noise floor only (chosen)** | 356 | 25 | 60 | 194 µm² |
| 54 µm³ (0.2 of expected) | 313 | 13 | 53 | **263 µm²** |

The spread halves, but 7 of 60 synapses disappear and the largest contact *grows*. Checking the
removed objects one at a time explains it — they are not debris:

| object | volume | overlap in x,y with a full nucleus | gap in z |
|---|---|---|---|
| 3 | 15 µm³ | 0.99 | 9 planes |
| 11 | 50 µm³ | 1.00 | 7 planes |
| 17 | 45 µm³ | 1.00 | 5 planes |
| 4, 9, 14, 15 | 12–51 µm³ | 0.00 | — |

A piece of a nucleus that the z-stitching failed to join would overlap in x and y *and* sit in an
adjacent plane. Nothing does. The objects with perfect x,y overlap are 5–9 planes away — they are
cells sitting above or below another cell. The rest do not overlap anything; 14 and 15 are 46 and
51 µm³ against a 54 µm³ cut-off, i.e. ordinary nuclei that Cellpose drew slightly small.

A cell whose only seed is removed has no marker left, so the watershed gives its whole body to
the neighbour and their shared border is reported as one large false contact. Hence
`--min-nucleus-fraction` stays at a floor against single specks (0.075 ≈ 20 µm³ for an 8 µm
nucleus) rather than becoming a real filter.

### Merging small objects into the nucleus beside them

Follows from the above: rather than dropping an undersized object, merge it into a full nucleus
within 4 µm, so pieces rejoin their parent while isolated small cells keep their seed. It does
not work — in a packed field every small object is within 4 µm of *some* nucleus, including ones
belonging to other cells. The problem pair stayed at 249 µm² against 194 µm² without it. Removed
again; the measurement above shows there were no pieces to rejoin in the first place.

### A larger Cellpose diameter, and two diameters

The two cell types have genuinely different nuclei — tumour cells median 9.2 µm across, T cells
6.3 µm — and Cellpose takes a single expected diameter. At 8 µm a tumour nucleus is sometimes cut
into three, and two nuclei are sometimes returned as one 14 µm object whose body then covers both
cells. Both failures are visible in `segmentation_qc/`.

A sweep of the diameter looked promising:

| diameter | tumour-side nuclei | median diameter | T-side nuclei |
|---|---|---|---|
| 7 µm | 5 | 6.3 µm | 24 |
| **8 µm (kept)** | 4 | 9.2 µm | 23 |
| 10 µm | 3 | 12.1 µm | 23 |
| 12 µm | 4 | 9.2 µm | 17 |
| 14 µm | 2 | 13.8 µm | 13 |

10 µm appears to hold the tumour nuclei together without costing T cells. Run end to end it does
not: fewer tumour objects meant *merged*, not *whole*, and the largest contact in that field went
from 225 to 386 µm². Counting objects is not a measure of correctness.

Detecting at two diameters and keeping each where it belongs fails the same way:

| attempt | largest contact in that field |
|---|---|
| single diameter 8 µm (kept) | 225 µm² |
| two diameters, chosen per object by GFP | 387 µm² |
| two diameters, each clipped to its side of the GFP boundary | 349 µm² |
| single diameter 10 µm | 386 µm² |

At 12 µm a tumour detection swallows the T cell next to it, and clipping at the GFP boundary
then fragments the T cell nuclei instead (bodies of 191, 368, 376 µm³ where a T cell is ~1200).

### Forcing GFP+ voxels out of tumour bodies

The all-channel footprint removed the old guarantee that a GFP+ voxel belongs to a T cell, so the
obvious repair is to hand every GFP+ voxel in a tumour label to the nearest T cell. The T cells
then inflate — one reached 2051 µm³ against a median of ~1200 — and the contact grew with them,
to 386 µm² in the field used for testing. Rejected.

The common thread in these four: **once two nuclei have been detected as one, the information
that there were two cells is gone, and no later step recovers it.**

### Flagging the fused nuclei instead of fixing them

If a fused nucleus cannot be split, the next best thing is to mark the rows it affects. Two
criteria were tried, neither of which works.

| criterion | nuclei flagged | synapses flagged | verdict |
|---|---|---|---|
| volume > 2x the median for that cell type | T 13 %, tumour 27 % | 30 of 60 (50 %) | useless — every cell here is clipped in z, so most of the volume spread is truncation, not fusion |
| x,y footprint > 2x the median (insensitive to z truncation) | T 1 %, tumour 8 % | 3 of 60 (5 %) | plausible rate, but it misses the cells that are actually wrong: none of the five largest contacts is flagged, and neither are the two nuclei in `4_stack5` that are visibly a Raji cut into three and a Raji fused with a T cell |

A flag that misses the obvious failures is worse than none, because an unflagged row then reads
as checked. Both were removed. The fused nuclei stay in the table unmarked, and
`segmentation_qc/` remains the only way to find them — look for two cells inside one outline, or
a line drawn through one cell.

---

## Why there is a manual step

Five attempts above failed for the same reason: the seeds are wrong, and nothing downstream can
tell that two cells were detected as one. The bodies grown from *correct* seeds are good, so the
cheapest reliable fix is to let a person place the seeds — one point per cell and its type — and
leave everything else automatic. That is `tools/curate_seeds.ijm` and `--curation`, described in
the README. It is not a workaround for a bug in the code; it is where the information the images
no longer carry has to come from.

## What would actually fix this

Not a parameter. The DAPI channel has to be smoothed by 1 µm before the nuclei can be detected at
all, which means the boundary between two touching nuclei is not in the image to begin with. A
brighter DAPI exposure, and a finer z step so a nucleus is sampled over more than 4–5 planes,
would remove the fusion and the fragmentation together — and would also make nucleus size
trustworthy enough for the debris filter above to become usable. See **Next experiments** in the
[README](../README.md).
