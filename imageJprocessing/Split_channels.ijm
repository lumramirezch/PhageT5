//Macro Split_channels.ijm
// Copy the directory tree with the following command on a linux-based system:
// rsync -a -f"+ */" -f"- *" path/to/source/ path/to/destination/
/*
 * This script will help you to split recursively czi files  
 * as a single tiff file per channel.
 */
//Path to the input and output
wfolder1 = "C:/004-Results/004-Microscopy/2020070301_03-OR63-T5-timecourse/" ; 
wfolder2I = "Data" ;
wfolder2O = "split" ;
wfolder3 = "/202007" ;
gene="03-T30-05-amA2" ;

// This is the input, images in czi 
directory = wfolder1 + wfolder2I + wfolder3 + gene + "/";

// Folder for the output
outputDirectory = wfolder1 + wfolder2O + wfolder3 + gene+"/";

// Declare variables
oneFilePerSlice = true;
fileList = getFileList(directory);

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

