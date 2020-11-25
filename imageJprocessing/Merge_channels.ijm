//Merge tiff files
// To get the best out of this script, copy the directory tree with 
// the following command on a linux-based system:
// rsync -a -f"+ */" -f"- *" path/to/source/ path/to/destination/
/*
 * This script will help you to crop 1/25 of images at a specific position, 
 * merge the images, add a 20-px scale bar, flatten
 * the merged and single images to save as tiff.
 */
//Path to the input and output
wfolder1 = "C:/004-Results/004-Microscopy/2020070301_03-OR63-T5-timecourse/" ; 
wfolder2I = "split" ;
wfolder2O = "merged" ;
wfolder3 = "/202007" ;
gene="03-T30-05-amA2" ;

// x position to crop the image
roiX=2*zoom; 
// y position to crop the image
roiY=2*zoom;


// Input directory (split tiff images)
outputDirectory = wfolder1 + wfolder2I + wfolder3 + gene + "/";

// Output directory (merged and cropped images)
outputMerged = wfolder1 + wfolder2O + wfolder3 + gene + "/";

//declare variables
input = true;
n=true;

//Zoom
zoom=1/5;

//Crop, add a white scale bar, flatten, and save the image
function cropFileW(input,n){
	open(outputDirectory + input);
	pixelH = getHeight();
	pixelW = getWidth();
	makeRectangle(pixelW*roiX, pixelH*roiY, pixelH*zoom, pixelW*zoom);
	run("Crop"); 
	setColor(255, 255, 255);
	fillRect(getHeight()/(1+zoom),getWidth()/(1+zoom),20,5);
	run("Enhance Contrast", "saturated=0.20 normalize");
	saveAs("tiff", outputMerged + substring(input,0,input.length-7)+"-"+n);
	close("*");
	}

//Crop and put a black scale bar
function cropFileB(input,n){
	open(outputDirectory + input);
	pixelH = getHeight();
	pixelW = getWidth();
	makeRectangle(pixelW*roiX, pixelH*roiY, pixelH*zoom, pixelW*zoom);
	run("Crop");
	//run("Enhance Contrast", "saturated=0.20 normalize");
	setColor(0, 0, 0);
	fillRect(getHeight()/(1+zoom),getWidth()/(1+zoom),20,5);
	saveAs("tiff", outputMerged + substring(input,0,input.length-7)+"-"+n);
	close("*");
	}
//merge two channels
function mergeFile(input2, input3) {
	open(outputDirectory + input2);
	open(outputDirectory + input3);
	// c1 is for Red and c3 is for blue
	run("Merge Channels...", "c3=[" + input2 + "] c1=[" + input3 + "] keep ignore");
	pixelH = getHeight();
	pixelW = getWidth();
	makeRectangle(pixelW*roiX, pixelH*roiY, pixelH*zoom, pixelW*zoom);
	run("Crop");
	setColor(255, 255, 255);
	run("Enhance Contrast", "saturated=0.20 normalize");
	fillRect(getHeight()/(1+zoom),getWidth()/(1+zoom),20,5);
	saveAs("tiff", outputMerged + substring(input2,0,input2.length-7));
	close("*");
	}


setBatchMode(true);

fileList = getFileList(outputDirectory);
fileList2 = List.clear();


input0 = fileList[1];
input1 = fileList[2];
input2 = fileList[3];

for (i=0; i<fileList.length; i++) {
	// ==0 is for phase contrast, ==1 for FM4-64 (red), and ==2 for DAPI (blue)
	if (substring(fileList[i],fileList[i].length-6,fileList[i].length-5)==0){
		input0 = fileList[i];
		cropFileW(input0,0);
		};
	if (substring(fileList[i],fileList[i].length-6,fileList[i].length-5)==1){
		input1 = fileList[i];
		cropFileW(input1,1);
		};
	if (substring(fileList[i],fileList[i].length-6,fileList[i].length-5)==2){
		input2=fileList[i];
		cropFileB(input2,2);
		};
	// if the name of both files is the same, they will be merged
	if (substring(input1,0,input1.length-7)==substring(input2,0,input2.length-7)){
		mergeFile(input0, input1);
		};
	}

showStatus("Finished.");
setBatchMode(false);