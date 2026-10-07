// Fiji/ImageJ macro: correct the cells the pipeline found, by hand, for a whole folder.
//
// The automatic detection gets the cell *bodies* right once it has the right seeds, and gets
// the seeds wrong in three ways: one nucleus found twice, two nuclei found as one, or the wrong
// cell type. All three are fixed here, by editing one marker per cell.
//
//   Fiji > Plugins > Macros > Run...  and pick this file
//   1st dialog = the folder of images that were analysed (.czi or *_stack.tif)
//   2nd dialog = the results folder of the automatic run (it holds seeds/ and rois/)
//   3rd dialog = where to save the corrected markers (the repository's curation/ folder)
//   then a short options dialog
//
// It walks through the images one after another, two rounds each, T cells then tumour cells.
// In every round the cells are both circles on the image and rows in the ROI Manager, the
// window listing them beside the image. Each circle carries its row number.
//
//   to ADD a cell     click it on the image, then press  t
//   to REMOVE a cell  click its row in the ROI Manager and press Delete
//                     (ctrl-click or shift-click to take several rows at once)
//   to MOVE a cell    drag its circle
//
// The number drawn on a circle is the name of its row, so the two always point at the same
// cell - including after deletions, which leave gaps in the numbering rather than shifting it.
//
// The plane does not matter: click anywhere on the cell, on whichever plane you like. Each
// marker is placed at the plane where its own DAPI column is brightest.
//
// Each image is saved as soon as its two rounds are done, so stopping half way loses nothing:
// leave "skip images already corrected" ticked next time and carry on where you left off.
//
// Then re-run the analysis with  --curation <the folder you saved to>.

var editedY = newArray(0);   // second return value of editRound()
var markerPx = 26;           // diameter of the circles, set per image from the pixel size

imageDir = getDirectory("The folder of images that were analysed");
resultsDir = getDirectory("The results folder of the automatic run");
outDir = getDirectory("Where to save the corrected markers (curation/)");

Dialog.create("Correct the cells");
Dialog.addString("Only file names containing", "_stack", 20);
Dialog.addNumber("Marker size (um)", 5);
Dialog.addCheckbox("Skip images already corrected", true);
Dialog.addCheckbox("Skip images with no automatic result", true);
Dialog.show();
filter = Dialog.getString();
markerUm = Dialog.getNumber();
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

run("ROI Manager...");
run("Labels...", "color=white font=18 show use draw");

for (n = 0; n < todo.length; n++) {
    name = stripExtension(todo[n]);
    where = "(" + (n + 1) + "/" + todo.length + ") " + name;
    openImage(imageDir + todo[n]);
    id = getImageID();
    getPixelSize(unit, pw, ph);
    if (startsWith(unit, "micro") || unit == "um") markerPx = round(markerUm / pw);
    else markerPx = 26;
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
    roiManager("reset");
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

// the automatic cell outlines, as a passive overlay so the ROI Manager stays free for the cells
function showOutlines(roiPath) {
    if (!File.exists(roiPath)) return;
    roiManager("reset");
    roiManager("Open", roiPath);
    roiManager("Show All without labels");
    run("From ROI Manager");
    roiManager("reset");
    Overlay.setStrokeColor("white");
}

// One round: every cell of one type is a numbered circle on the image and a row in the ROI
// Manager. Returns the edited x, and sets editedY to the matching y.
function editRound(id, xs, ys, where, what, colour) {
    selectImage(id);
    run("Select None");
    roiManager("reset");
    for (i = 0; i < xs.length; i++) {
        makeOval(xs[i] - markerPx / 2, ys[i] - markerPx / 2, markerPx, markerPx);
        roiManager("add");
    }
    if (roiManager("count") > 0) {
        roiManager("Deselect");
        roiManager("Set Color", colour);
        roiManager("Set Line Width", 2);
        // name every cell after its number, so the row in the list and the number drawn on the
        // image are the same text. ImageJ's own names are slice-and-coordinate strings, which
        // match nothing on screen, and its index labels renumber themselves after a deletion.
        for (i = 0; i < roiManager("count"); i++) {
            roiManager("select", i);
            roiManager("rename", "" + (i + 1));
        }
        roiManager("Deselect");
    }
    run("Labels...", "color=white font=18 show use draw");
    roiManager("Show All with labels");
    run("Select None");
    setTool("point");
    run("Point Tool...", "type=Hybrid color=" + colour + " size=Large");

    waitForUser(where + " - the " + what,
        "ADD     click the cell, then press  t\n"
        + "REMOVE  click its row in the ROI Manager, press Delete\n"
        + "MOVE    drag the circle\n \n"
        + "The plane does not matter.\n \n"
        + "OK = done with the " + what + ".   Cancel = stop (finished images are saved).");

    // a cell clicked but not yet committed with 't' - take it rather than lose it
    if (selectionType() == 10 && roiManager("index") == -1)
        roiManager("add");

    keptX = newArray(0); keptY = newArray(0);
    for (i = 0; i < roiManager("count"); i++) {
        roiManager("select", i);
        if (selectionType() == 10) {           // a point, or several added in one go
            getSelectionCoordinates(px, py);
            for (k = 0; k < px.length; k++) {
                keptX = Array.concat(keptX, px[k]);
                keptY = Array.concat(keptY, py[k]);
            }
        } else {                                // a circle: its centre is the cell
            Roi.getBounds(bx, by, bw, bh);
            keptX = Array.concat(keptX, round(bx + bw / 2));
            keptY = Array.concat(keptY, round(by + bh / 2));
        }
    }
    roiManager("reset");
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
