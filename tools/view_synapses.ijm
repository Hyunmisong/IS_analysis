// Fiji/ImageJ macro: look at what the pipeline found, on the original image.
//
// Loads the cell outlines and the detected synapses as ROIs on top of the analysed stack:
//   cyan   = T cell body        (name z07_T3  = plane 7, T cell 3)
//   yellow = tumour cell body   (name z07_Tumour12)
//   white  = the synapse itself (name z07_IS01 = the contact measured as synapse 1)
// The numbers are the T_Cell / Tumour_Cell / IS columns of results/synapse_table.csv.
//
//   Fiji > Plugins > Macros > Run...  and pick this file
//   first dialog  = the analysed image (*_stack.tif, or the original .czi)
//   second dialog = the results folder the pipeline wrote
//
// In the ROI Manager, untick "Show All" to clear the image, or select a single ROI to see one
// cell. Scrolling through z shows only the ROIs of the current plane.

stackPath = File.openDialog("The analysed image (*_stack.tif or .czi)");
resultsDir = getDirectory("The results folder");
name = stripExtension(File.getName(stackPath));
roiPath = resultsDir + "rois" + File.separator + name + ".zip";

if (!File.exists(roiPath))
    exit("No ROI file for this image:\n" + roiPath + "\n\nRun the analysis on it first.");

if (endsWith(toLowerCase(stackPath), ".czi"))
    run("Bio-Formats Importer",
        "open=[" + stackPath + "] color_mode=Composite view=Hyperstack stack_order=XYCZT");
else
    open(stackPath);
Stack.setDisplayMode("composite");
run("Enhance Contrast", "saturated=0.35");

roiManager("reset");
roiManager("Open", roiPath);
roiManager("Show All with labels");
print(name + ": " + roiManager("count") + " outlines loaded from " + roiPath);

function stripExtension(fileName) {
    dot = lastIndexOf(fileName, ".");
    if (dot < 0) return fileName;
    return substring(fileName, 0, dot);
}
