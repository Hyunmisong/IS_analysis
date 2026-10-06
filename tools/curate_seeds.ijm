// Fiji/ImageJ macro: correct the cells the pipeline found, by hand, for a whole folder.
//
// The automatic detection gets the cell *bodies* right once it has the right seeds, and gets
// the seeds wrong in three ways: one nucleus found twice, two nuclei found as one, or the wrong
// cell type. All three are fixed here, by editing one point per cell.
//
//   Fiji > Plugins > Macros > Run...  and pick this file
//   1st dialog = the folder of images that were analysed (.czi or *_stack.tif)
//   2nd dialog = the results folder of the automatic run (it holds seeds/ and rois/)
//   3rd dialog = where to save the corrected points (the repository's curation/ folder)
//   then a short options dialog
//
// It then walks through the images one after another. For each one you get two rounds, T cells
// and tumour cells. In each round:
//   * the points already found are shown, and the automatic cell outlines sit underneath as a
//     white overlay so you can see what the program thought;
//   * the multi-point tool is selected. Click to add a cell. Alt-click exactly on a point to
//     remove it - that only works within a few pixels, so zoom in (+) if you are missing;
//   * then a second step asks you to draw around any points that should go. That is the
//     reliable way to delete: draw a rectangle or a freehand loop round them and click OK.
//     Leave the image unselected to keep everything;
//   * to change a cell's type, remove it in one round and add it in the other;
//   * z does not matter - click anywhere on the cell, on whichever plane you like. The analysis
//     puts each point on the plane where its own DAPI signal is brightest;
//   * click OK to go on to the next round, or Cancel to stop for the day.
//
// Each image is saved as soon as its two rounds are done, so stopping half way loses nothing:
// leave "skip images already corrected" ticked next time and carry on where you left off.
//
// Then re-run the analysis with  --curation <the folder you saved to>.

var editedY = newArray(0);   // second return value of editRound()

imageDir = getDirectory("The folder of images that were analysed");
resultsDir = getDirectory("The results folder of the automatic run");
outDir = getDirectory("Where to save the corrected points (curation/)");

Dialog.create("Correct the cells");
Dialog.addString("Only file names containing", "_stack", 20);
Dialog.addCheckbox("Skip images already corrected", true);
Dialog.addCheckbox("Skip images with no automatic result", true);
Dialog.show();
filter = Dialog.getString();
skipDone = Dialog.getCheckbox();
skipMissing = Dialog.getCheckbox();

list = getFileList(imageDir);
todo = newArray(0);
for (i = 0; i < list.length; i++) {
    if (!isImage(list[i]) || indexOf(list[i], filter) < 0) continue;
    name = stripExtension(list[i]);
    if (skipDone && File.exists(outDir + name + ".csv")) continue;
    if (skipMissing && !File.exists(resultsDir + "seeds" + File.separator + name + ".csv"))
        continue;
    todo = Array.concat(todo, list[i]);
}
if (todo.length == 0)
    exit("Nothing to do in " + imageDir + "\n(filter \"" + filter + "\", and the skip options).");

for (n = 0; n < todo.length; n++) {
    name = stripExtension(todo[n]);
    where = "(" + (n + 1) + "/" + todo.length + ") " + name;
    openImage(imageDir + todo[n]);
    id = getImageID();
    showOutlines(resultsDir + "rois" + File.separator + name + ".zip");

    tx = newArray(0); ty = newArray(0); ux = newArray(0); uy = newArray(0);
    seedPath = resultsDir + "seeds" + File.separator + name + ".csv";
    if (File.exists(seedPath)) {
        lines = split(File.openAsString(seedPath), "\n");
        for (i = 1; i < lines.length; i++) {           // row 0 is the header
            row = split(replace(lines[i], "\r", ""), ",");
            if (row.length < 3) continue;
            kind = replace(toLowerCase(row[2]), " ", "");
            if (startsWith(kind, "tu")) {
                ux = Array.concat(ux, parseInt(row[0])); uy = Array.concat(uy, parseInt(row[1]));
            } else {
                tx = Array.concat(tx, parseInt(row[0])); ty = Array.concat(ty, parseInt(row[1]));
            }
        }
    }

    tx = editRound(id, tx, ty, where, "T cells (GFP+)", "cyan");
    ty = editedY;
    ux = editRound(id, ux, uy, where, "tumour cells (GFP-)", "yellow");
    uy = editedY;

    out = "x,y,type\n";
    for (i = 0; i < tx.length; i++) out = out + tx[i] + "," + ty[i] + ",T\n";
    for (i = 0; i < ux.length; i++) out = out + ux[i] + "," + uy[i] + ",Tumour\n";
    File.saveString(out, outDir + name + ".csv");
    print(where + ": " + tx.length + " T cells, " + ux.length + " tumour cells saved");
    selectImage(id);
    close();
}
print("Done. Now re-run:");
print("  python scripts/run_analysis.py ... --curation \"" + outDir + "\" --out results_curated");

function openImage(path) {
    if (endsWith(toLowerCase(path), ".czi"))
        run("Bio-Formats Importer",
            "open=[" + path + "] color_mode=Composite view=Hyperstack stack_order=XYCZT");
    else
        open(path);
    Stack.setDisplayMode("composite");
    run("Enhance Contrast", "saturated=0.35");
}

// the automatic outlines, as a passive overlay so the ROI Manager stays free for the points
function showOutlines(roiPath) {
    if (!File.exists(roiPath)) return;
    roiManager("reset");
    roiManager("Open", roiPath);
    roiManager("Show All without labels");
    run("From ROI Manager");
    roiManager("reset");
    Overlay.setStrokeColor("white");
}

// Show one cell type's points, wait for the user to edit them, return the edited x (and set
// editedY). ImageJ macros return one value, hence the global for the second array.
// Show one cell type's points, let the user edit them, return the edited x (and set editedY).
// ImageJ macros return one value, hence the global for the second array.
function editRound(id, xs, ys, where, what, colour) {
    selectImage(id);
    run("Select None");
    setTool("multipoint");
    // Hybrid = a cross with a dot in the middle, which stays visible over a noisy image, and the
    // largest size also widens the few-pixel target that alt-click has to hit.
    run("Point Tool...", "type=Hybrid color=" + colour + " size=[Extra Large] label");
    if (xs.length > 0) makeSelection("point", xs, ys);
    waitForUser(where + " - mark the " + what,
        "Click to add a cell.\n"
        + "Alt-click exactly on a point to remove it (zoom in with + if you keep missing).\n"
        + "The plane does not matter.\n\n"
        + "OK = go on to deleting.   Cancel = stop (finished images are saved).");
    if (selectionType() == 10)
        getSelectionCoordinates(xs, ys);
    else
        { xs = newArray(0); ys = newArray(0); }

    xs = deleteRound(id, xs, ys, where, what, colour);
    ys = editedY;
    editedY = ys;
    return xs;
}

// The reliable way to remove points: draw a region around them. Alt-click has to land within a
// few pixels of a point, which is hard on a crowded field; a loop around them never misses.
function deleteRound(id, xs, ys, where, what, colour) {
    selectImage(id);
    run("Select None");
    marks = Overlay.size;
    for (i = 0; i < xs.length; i++) {
        makePoint(xs[i], ys[i], "small " + colour + " hybrid");
        Overlay.addSelection;
    }
    run("Select None");
    setTool("rectangle");
    waitForUser(where + " - remove any wrong " + what,
        "Draw a rectangle or a freehand loop around the points to delete.\n"
        + "Every point inside it goes.\n\n"
        + "OK with nothing drawn = keep them all.");
    keptX = newArray(0); keptY = newArray(0);
    area = (selectionType() >= 0 && selectionType() <= 4) || selectionType() == 9;
    for (i = 0; i < xs.length; i++) {
        if (area && selectionContains(xs[i], ys[i])) continue;
        keptX = Array.concat(keptX, xs[i]);
        keptY = Array.concat(keptY, ys[i]);
    }
    for (i = Overlay.size - 1; i >= marks; i--)   // take the temporary marks back off
        Overlay.removeSelection(i);
    run("Select None");
    editedY = keptY;
    return keptX;
}

function isImage(fileName) {
    lower = toLowerCase(fileName);
    return endsWith(lower, ".czi") || endsWith(lower, ".tif") || endsWith(lower, ".tiff");
}

function stripExtension(fileName) {
    dot = lastIndexOf(fileName, ".");
    if (dot < 0) return fileName;
    return substring(fileName, 0, dot);
}
