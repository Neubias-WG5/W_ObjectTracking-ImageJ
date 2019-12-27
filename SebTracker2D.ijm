// Author: Sébastien Tosi (IRB Barcelona)
// Contact: sebastien.tosi@irbbarcelona.org
// Version: 1.0
// Date: 17/12/2019

// Path to input image and results
inputDir = "C:\\Users\\stosi\\Desktop\\in";
outputDir = "C:\\Users\\stosi\\Desktop\\out";

// Parameters for CTC time-lapse
medrad = 5;		  // Median filter radius (pix)
thr= 105;		  // Global threshold (fixed)
erodrad = 5;		  // Erosion radius to seed objects (pix)
dmapds = 2;    		  // Used to speed up (set from 1-3, recommended 2)
noisetol = 3;		  // Noise tolerance to split out touching objects (distance map gray levels)

// Read arguments from command line
arg = getArgument();
parts = split(arg, ",");
for(i=0; i<parts.length; i++) 
{
	nameAndValue = split(parts[i], "=");
	if (indexOf(nameAndValue[0], "input")>-1) inputDir=nameAndValue[1];
	if (indexOf(nameAndValue[0], "output")>-1) outputDir=nameAndValue[1];
	if (indexOf(nameAndValue[0], "medrad")>-1) medrad=nameAndValue[1];
	if (indexOf(nameAndValue[0], "thr")>-1) thr=nameAndValue[1];
	if (indexOf(nameAndValue[0], "erodrad")>-1) erodrad=nameAndValue[1];
	if (indexOf(nameAndValue[0], "dmapds")>-1) dmapds=nameAndValue[1];
	if (indexOf(nameAndValue[0], "noisetol")>-1) noisetol=nameAndValue[1];
}

images = getFileList(inputDir);
for(img=0; img<images.length; img++) 
{

setBatchMode(true);

// Division file (should be stored to some temporary folder accessible for metrics computation)
FileName = images[img];
FirstDot = indexOf(FileName,".");
DivFile = outputDir+File.separator+substring(FileName,0,FirstDot)+".txt";

// Open image
open(inputDir+File.separator+FileName);

// Initialization
if(isOpen("ROI Manager"))
{
	selectWindow("ROI Manager");
	run("Close");
}
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
run("Set Measurements...", "area mean centroid stack redirect=None decimal=3");
run("Select None");

// Pre-processing
if(medrad>0)run("Median...", "radius="+d2s(medrad,0)+" stack");
setThreshold(thr, 65535);
setOption("BlackBackground", false);
run("Convert to Mask", "method=Default background=Dark");
run("Fill Holes", "stack");
if(noisetol>-1)
{
	if(noisetol==0)run("Watershed", "stack");
	else modWaterhsed(noisetol,dmapds);
}
if(erodrad>0)run("Minimum...", "radius="+d2s(erodrad,0)+" stack");
rename("SeedMask.tif");

// Object analysis (Step 1: find connection tubes)
run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
rename("Labels.tif");selectImage("SeedMask.tif");close();selectImage("Labels.tif");

// Object analysis (Step 2: Cut tubes at division points and reanalyze connection tubes)
run("Duplicate...", "title=SplitMask.tif duplicate");
Stack.getStatistics(voxelCount, mean, min, NObjs);
DifCounts = newArray(NObjs+1);
LastCounts = newArray(NObjs+1);
for(i=nSlices;i>=1;i--)
{
	ActiveCounts = newArray(NObjs+1);
	selectImage("SplitMask.tif");
	run("Select None");
	setSlice(i);
	wait(10);
	setThreshold(1,65535);
	run("Analyze Particles...", "  display clear slice");
	resetThreshold();
	for(j=0;j<nResults;j++)ActiveCounts[getResult("Mean",j)] = ActiveCounts[getResult("Mean",j)]+1;
	if(i<nSlices)for(j=1;j<NObjs+1;j++)DifCounts[j] = LastCounts[j] - ActiveCounts[j]; 
	for(j=1;j<NObjs+1;j++)
	{
		if(DifCounts[j]>0)
		{
			for(k=0;k<nResults;k++)
			{
				if(getResult("Mean",k)==j)
				{
					doWand(getResult("X",k),getResult("Y",k));
					roiManager("add");
					run("Set...", "value=0 slice");
				}	
			}
		}	
	}
	LastCounts = ActiveCounts;
}
selectImage("SplitMask.tif");
setThreshold(1,65535);
run("Convert to Mask", "method=Default background=Dark");
run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
rename("Labels-final.tif");
Stack.getStatistics(voxelCount, mean, min, NObjs2);
run("Select None");

// Restore connection tubes cuts
selectImage("Labels-final.tif");
for(i=0;i<roiManager("count");i++)
{
	roiManager("select",i);
	Slice = getSliceNumber();
	if(Slice>1)
	{
		setSlice(Slice-1);
		getStatistics(voxelCount, mean, min, ParentLbl);
		setSlice(Slice);
		run("Set...", "value="+d2s(ParentLbl,0)+" slice");
	}
}

// Analyze objects life span
run("Select None");
StartFrame = newArray(NObjs2+1);
EndFrame = newArray(NObjs2+1);
for(j=0;j<=NObjs2;j++)StartFrame[j] = -1;
for(i=0;i<nSlices;i++)
{
	setSlice(i+1);
	getHistogram(values, counts, 65536);
	for(j=1;j<=NObjs2;j++)
	{
		if(counts[j]>0)
		{
			if(StartFrame[j]==-1)StartFrame[j] = i;
			EndFrame[j] = i;
		}
	}
}
setSlice(1);

// Create division file (only objects with two largest areas are considered children, other object is appearing)
Str2 = "";
Marked = newArray(NObjs2+1);
for(i=0;i<roiManager("count");i++)
{	
	roiManager("select",i);
	Slice = getSliceNumber();
	getStatistics(voxelCount, mean, min, ParentLbl);
	setSlice(Slice+1);
	getHistogram(values, counts, 65536);
	Child1Lbl = 0;Cnt1 = 0;
	Child2Lbl = 0;Cnt2 = 0;
	for(j=1;j<65536;j++)
	{
		if(counts[j]>Cnt1)
		{
			if(Cnt1>0)
			{
				Child2Lbl = Child1Lbl;
				Cnt2 = Cnt1;
			}
			Child1Lbl = j;
			Cnt1 = counts[j];
		}
		else
		{
			if(counts[j]>Cnt2)
			{
				Child2Lbl = j;
				Cnt2 = counts[j];
			}
		}
	}
	//print("In frame "+d2s(getSliceNumber(),0)+" object "+d2s(ParentLbl,0)+" divides into object "+d2s(Child1Lbl,0)+" and object "+d2s(Child2Lbl,0));
	// Label / Start slice / End slice / ParentLbl
	Str2 = Str2 + d2s(Child1Lbl,0)+" "+d2s(StartFrame[Child1Lbl],0)+" "+d2s(EndFrame[Child1Lbl],0)+" "+d2s(ParentLbl,0)+"\n";
	Str2 = Str2 + d2s(Child2Lbl,0)+" "+d2s(StartFrame[Child2Lbl],0)+" "+d2s(EndFrame[Child2Lbl],0)+" "+d2s(ParentLbl,0)+"\n";
	Marked[Child1Lbl] = 1;
	Marked[Child2Lbl] = 1;
}
// Add dividing objects without parents to division file
Str = "";
for(j=1;j<=NObjs2;j++)if(Marked[j]==0)Str = Str + d2s(j,0)+" "+d2s(StartFrame[j],0)+" "+d2s(EndFrame[j],0)+" 0\n";

// Dilate to restore original objects shape
run("Select None");
if(erodrad>0)run("Maximum...", "radius="+d2s(erodrad,0)+" stack");

// Apply LUT
run("Rainbow RGB");
selectImage("Labels-final.tif");
setSlice(nSlices);resetMinAndMax();setSlice(1);

// Save label mask
//run("Bio-Formats Exporter", "save="+outputDir+File.separator+FileName+" compression=Uncompressed");
save(outputDir+File.separator+FileName);

// Save division file
File.saveString(Str+Str2, DivFile);

// Cleanup
run("Close All");
setBatchMode("exit & display");

}


// Modified watershed (find significant regional maxima)
function modWaterhsed(noisetol,dmapds)
{
	rename("Mask.tif");
	W = getWidth;H = getHeight();D = nSlices;
	run("Scale...", "x="+d2s(1/dmapds,2)+" y="+d2s(1/dmapds,2)+" z=1.0 interpolation=None process create title=Dmap");
	run("Distance Map", "stack");
	N = nSlices;
	for(i=0;i<N;i++)
	{
		selectImage("Dmap");
		setSlice(i+1);
		run("Find Maxima...", "noise="+d2s(round(noisetol/dmapds),1)+" output=[Segmented Particles] light");
		if(i==0)rename("WatershedLinesDS");
		else 
		{
			rename("WatershedLinesDS_");
			run("Concatenate...", "  title=WatershedLinesDS image1=WatershedLinesDS image2=WatershedLinesDS_ image3=[-- None --]");
		}
		showProgress(i/N);
	}
	run("Scale...", "x="+d2s(dmapds,2)+" y="+d2s(dmapds,2)+" z=1 width="+d2s(W,0)+" height="+d2s(H,0)+" depth="+d2s(D,0)+" interpolation=None process create title=WatershedLines");
	selectImage("Dmap");
	close();
	selectImage("WatershedLinesDS");
	close();
	imageCalculator("AND stack", "Mask.tif","WatershedLines");
	selectImage("WatershedLines");
	close();
}
