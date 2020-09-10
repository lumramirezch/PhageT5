//Macro Split_channels.ijm
// Number and Name of the gene
gene="18-A1A2"

// This is the input, images in czi 
directory = "./2020090301-OR63-ectopic-140min/Data/20200903-ectopic-140min-"+gene+"/";

// Folder for the output
outputDirectory = "./2020090301-OR63-ectopic-140min/split/20200903-ectopic-140min-"+gene+"/";

fileList = getFileList(directory);

// directory = getDirectory(input);
fileList = getFileList(directory);

// outputDirectory = getDirectory(output);

run("Bio-Formats Macro Extensions");
setBatchMode(true);

for (i=0; i<fileList.length; i++) {
  file = directory + fileList[i];
  if (oneFilePerSlice) {
    Ext.setId(file);
    Ext.getImageCount(imageCount);

    for (image=0; image<imageCount; image++) {
      Ext.openImage("czi", image);
      outFile = outputDirectory + substring(fileList[i], 0, fileList[i].length-4) + "-" + image +".tiff";
      saveFile(outFile);
      close();
    }
    Ext.close();
  }
  else {
    Ext.openImagePlus(file);
    outFile = outputDirectory + substring(fileList[i], 0, fileList[i].length-4) + ".tiff";
    saveFile(outFile);
    close();
  }
}

showStatus("Finished.");
setBatchMode(false);

function saveFile(outFile) {
   run("Bio-Formats Exporter", "save=[" + outFile + "] compression=Uncompressed");
}

