// Fiji/ImageJ macro: inspect one result of the pipeline.
//
// Opens a *_stack.tif side by side with the cell-body label image the pipeline wrote to
// results/masks/, with the two windows synchronised, so the segmentation can be scrolled
// through plane by plane in Fiji.
//
//   Fiji > Plugins > Macros > Run...  and pick this file
//   first dialog  = the *_stack.tif that was analysed
//   second dialog = results/masks/<same name>_bodies.tif
//
// Every cell has its own colour in the label image, and its number is the T_Cell /
// Tumour_Cell column of results/synapse_table.csv.

stackPath = File.openDialog("The analysed *_stack.tif");
maskPath = File.openDialog("The matching results/masks/*_bodies.tif");

open(stackPath);
Stack.setDisplayMode("composite");
run("Enhance Contrast", "saturated=0.35");

open(maskPath);
run("glasbey on dark");
getStatistics(area, mean, min, max);
setMinAndMax(0, max);

run("Tile");
run("Synchronize Windows");
