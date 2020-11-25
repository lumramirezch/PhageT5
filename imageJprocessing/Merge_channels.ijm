//Merge tiff files

// Number and Name of the gene
gene="19-F";
// x position to crop the image
roiX=1.5/3; 
// y position to crop the image
roiY=1.8/3;

// Input to merge 
outputDirectory = "./2020090301-OR63-ectopic-140min/split/20200903-ectopic-140min-"+gene+"/";

// Folder for the output
outputMerged = "./2020090301-OR63-ectopic-140min/merged/20200903-ectopic-140min-"+gene+"/";

function cropFile(input1){
	open(outputDirectory + input1);
	pixelH = getHeight();
	pixelW = getWidth();
	makeRectangle(pixelW*roiX, pixelH*roiY, pixelH/3, pixelW/3);
	run("Crop");
	setColor(0, 0, 0);
	fillRect(getHeight()*5/6,getWidth()*5/6,50,5);
	run("Enhance Contrast", "saturated=0.20 normalize");
	saveAs("tiff", outputMerged + substring(input1,0,input1.length-7)+"-0");
	close("*");
	}

function mergeFile(input2, input3) {
	open(outputDirectory + input2);
	open(outputDirectory + input3);
	// c1 is for Red and c3 is for blue
	run("Merge Channels...", "c3=[" + input2 + "] c1=[" + input3 + "] keep ignore");
	pixelH = getHeight();
	pixelW = getWidth();
	makeRectangle(pixelW*roiX, pixelH*roiY, pixelH/3, pixelW/3);
	run("Crop");
	setColor(255, 255, 255);
	fillRect(getHeight()*5/6,getWidth()*5/6,50,5);
	run("Enhance Contrast", "saturated=0.20 normalize");
	saveAs("tiff", outputMerged + substring(input2,0,input2.length-7));
	close("*");
	}


setBatchMode(true);

fileList = getFileList(outputDirectory);
fileList2 = List.clear();


input1 = fileList[1];
input2 = fileList[2];
input3 = fileList[3];

for (i=0; i<fileList.length; i++) {
	// ==0 is for phase contrast, ==1 for FM4-64 (red), and ==2 for DAPI (blue)
	if (substring(fileList[i],fileList[i].length-6,fileList[i].length-5)==2){
		input1 = fileList[i];
		cropFile(input1);
		};
	if (substring(fileList[i],fileList[i].length-6,fileList[i].length-5)==0){
		input2 = fileList[i];
		};
	if (substring(fileList[i],fileList[i].length-6,fileList[i].length-5)==1){
		input3=fileList[i];
		};
	// if the name of both files is the same, they will be merged
	if (substring(input2,0,input2.length-7)==substring(input3,0,input3.length-7)){
		mergeFile(input2, input3);
		};
	}

showStatus("Finished.");
setBatchMode(false);