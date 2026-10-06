// Fiji/ImageJ macro: batch-convert confocal files (.czi, .lsm, .nd2, ...) into the
// *_stack.tif hyperstacks that the Python pipeline reads.
//
//   Fiji > Plugins > Macros > Run...  and pick this file
//   first dialog  = folder with the raw microscope files
//   second dialog = output folder (use data/raw of this repository)
//
// Channel order is kept as acquired. The pipeline expects DAPI, GFP, mCherry; if your
// acquisition order differs, pass it on the command line, e.g. --channels 2,1,0.

inDir = getDirectory("Folder with the raw microscope files");
outDir = getDirectory("Output folder (data/raw)");
exts = newArray(".czi", ".lsm", ".nd2", ".oib", ".oif", ".lif");

list = getFileList(inDir);
setBatchMode(true);
n = 0;
for (i = 0; i < list.length; i++) {
    if (!isImageFile(list[i], exts)) continue;
    run("Bio-Formats Importer",
        "open=[" + inDir + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
    name = stripExtension(list[i]);
    saveAs("Tiff", outDir + name + "_stack.tif");
    close();
    n++;
    print("saved " + name + "_stack.tif");
}
setBatchMode(false);
print("done: " + n + " file(s) -> " + outDir);

function isImageFile(name, exts) {
    lower = toLowerCase(name);
    for (k = 0; k < exts.length; k++)
        if (endsWith(lower, exts[k])) return true;
    return false;
}

function stripExtension(name) {
    dot = lastIndexOf(name, ".");
    if (dot < 0) return name;
    return substring(name, 0, dot);
}
