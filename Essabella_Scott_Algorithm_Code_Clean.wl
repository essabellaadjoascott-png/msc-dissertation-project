(* Input cell 1 *)
ClearAll["Global`*"]

(* Input cell 2 *)
(*Importing Images to Mathematica for Processing (corresponds to Section 3 (Data) in the thesis)
First, all of the images in the provided dataset were imported to Mathematica for processing (see below). The colour in each image name, for example "Blue" in sun2Blue, relates to the colour of the image. Specifically, Blue relates to the DAPI images, Red relates to the Alexa 594 Channel images, and Grey relates to the TRANS images (as described in section 3 of the thesis).*) 
 
 sun2 = Import["C:\\Users\\First_User\\OneDrive\\Desktop\\PROJECT\Images - proteins\\Sun2 stack.tif"];
sun2Blue=sun2[[1]]; 
 sun2Red =sun2[[2]]; 
 sun2Grey = sun2[[3]];
nesprin1 = Import["C:\\Users\\First_User\\OneDrive\\Desktop\\PROJECT\Images - proteins\\Nesprin 1 stack.tif"];
nesprin1Blue =nesprin1[[1]]; 
 nesprin1Red = nesprin1[[2]]; 
 nesprin1Grey = nesprin1[[3]]; 
laminB1 = Import["C:\\Users\\First_User\\OneDrive\\Desktop\\PROJECT\Images - proteins\\Lamin B1 stack.tif"]; 
 laminB1Blue=laminB1[[1]]; 
 laminB1Red=laminB1[[2]]; 
 laminB1Grey = laminB1[[3]]; 
 
 (* displaying the brightfield (TRANS) Images of all the proteins in a row*) 
 
 GraphicsRow[{sun2Grey, nesprin1Grey, laminB1Grey}, ImageSize->Large]

(* Input cell 3 *)
(*The Mean, Dimensions, Pixel Intensity Type and Channels of the grayscale images were checked*) 
 
 ImageMeasurements[#, "Mean"]&/@{sun2Grey, nesprin1Grey, laminB1Grey}
ImageDimensions[#]&/@{sun2Grey, nesprin1Grey, laminB1Grey}
ImageType[#]&/@{sun2Grey, nesprin1Grey, laminB1Grey}
ImageChannels[#]&/@{sun2Grey, nesprin1Grey, laminB1Grey}

(* Input cell 4 *)
(*Since ImageChannels confirmed that each image contained three channels (meaning they were recognised as RGB images), I converted them to Grayscale to ensure that single-channel intensity values could be calculated correctly*) 
 
 sun2Grayscale=ColorConvert[sun2Grey, "Grayscale"];
nesprin1Grayscale=ColorConvert[nesprin1Grey, "Grayscale"];
laminB1Grayscale=ColorConvert[laminB1Grey, "Grayscale"]; 
 
 ImageDimensions[#]&/@{sun2Grayscale, nesprin1Grayscale, laminB1Grayscale}
ImageChannels[#]&/@{sun2Grayscale, nesprin1Grayscale, laminB1Grayscale} 
 ImageType[#]&/@{sun2Grayscale, nesprin1Grayscale, laminB1Grayscale}

(* Input cell 5 *)

 (*Intensity histograms for each grayscale image above was applied using 256 bins to show the full range of possible pixel values. The histogram for all three slides showed a Gaussian-like (normal) distribution*) 
 
 GraphicsColumn[{ImageHistogram[sun2Grayscale, "Byte", FrameTicks->Automatic, PlotLabel->"Sun2"], ImageHistogram[nesprin1Grayscale, "Byte", FrameTicks->Automatic, PlotLabel->"Nesprin1"], ImageHistogram[laminB1Grayscale, "Byte", FrameTicks->Automatic, PlotLabel->"LaminB1"]}, ImageSize->Large] 
 
 (*Also the Mean and Standard deviation values of each grayscale image was calculated to show range or spread of intensity values. It was observed that the mean intensity & standard deviation of the SUN2 brightfield image was significantly higher than those of the Nesprin-1 and Lamin-B1 images *) 
 
 ImageMeasurements[#, "Mean"]&/@{sun2Grayscale, nesprin1Grayscale, laminB1Grayscale} 
 ImageMeasurements[#, "StandardDeviation"]&/@{sun2Grayscale, nesprin1Grayscale, laminB1Grayscale} 


(* Input cell 6 *)
(*Histogram Equalization was implemented using "HistogramTransform" function in Mathematica to transform the pixel values to enhance contrast and flatten the histogram of each slide*) 
 
 sun2GTransformed = HistogramTransform[sun2Grayscale]; 
 nesprin1GTransformed = HistogramTransform[nesprin1Grayscale]; 
 laminB1GTransformed = HistogramTransform[laminB1Grayscale]; 
 
 (*The transformed images & their histograms of each grayscale image were displayed. After histogram transformation, the means and standard deviations of all three slides were similar. See Section 3.2, which shows transformation & Histogram of SUN2 grayscale images*) 
 
 GraphicsRow[{sun2GTransformed, nesprin1GTransformed, laminB1GTransformed}, ImageSize->Large] 
 
 GraphicsColumn[{ImageHistogram[sun2GTransformed, "Byte", FrameTicks->Automatic, PlotLabel->"SUN2"], ImageHistogram[nesprin1GTransformed, "Byte", FrameTicks->Automatic], ImageHistogram[laminB1GTransformed, "Byte", FrameTicks->Automatic]}, ImageSize->Large]

(* Input cell 7 *)
(*However since it was difficult to ascertain which nuclei are intrafusal and extrafusal respectively, the inverse of the images were obtained to enhance contrast and to define cell boundaries. Inverse makes white pixels dark and causes dark pixels to be bright. All these correspond to Section 3.2*) 
 
 sun2GInverse=Image[1-sun2GTransformed];
nesprin1GInverse=Image[1-nesprin1GTransformed]; 
 laminB1GInverse=Image[1-laminB1GTransformed]; 

GraphicsRow[{sun2GInverse, nesprin1GInverse, laminB1GInverse}, ImageSize -> Large]

(* Input cell 8 *)
(*In an attempt to segment cells to help in spatial characterization of proteins, the original segmentation pipeline by Shelley Buchan was employed and applied to inverted images. 
 Output was bad as requirement for watershed algorithm was not met. Subsequently the segmentation algorithm was tweaked by skipping first functions that enhance contrast, gaussian smoothing over a large radius and exploiting MinDetect (for local minima) and MaxDetect (for local maxima) functions based on the intensity values of cell boundaries.
 See Supplementary Material Figure S1 to see failed attempt on segmentation of cell boundaries of SUN2 brightfield images*) 
 
 segmentation[img_]:=Module[{contrasted, brighter, meanFilter, binaryMask, masked, grayscale, gaussianFilter, maxValues, centroids, segmented}, contrasted=2-img; 
 brighter=ImageMultiply[img, contrasted]; 
 meanFilter=MeanFilter[brighter, 2]; 
 binaryMask=Binarize[meanFilter, FindThreshold[meanFilter, Method -> "Entropy"]]; 
 masked=ImageMultiply[binaryMask, meanFilter]; 
 grayscale=ColorConvert[masked, "Grayscale"]; 
 gaussianFilter=GaussianFilter[grayscale, 25]; 
 maxValues=MaxDetect[gaussianFilter, Padding -> 1]; 
 centroids=ComponentMeasurements[maxValues, "Centroid"]; 
 segmented=Colorize[WatershedComponents[GradientFilter[masked, 2], centroids[[All, 2]]]]] (*original Segmentation Algorithm developed by Shelley Buchan and designed for segmentation of nuclei from DAPI images*) 
 
 GraphicsRow[{segmentation[sun2Blue], segmentation[nesprin1Blue], segmentation[laminB1Blue]}] (*shows segmentation of nuclei. Section 4.2 corresponds to segmentation of nucle from DAPI images of SUN2 protein*)


(* Input cell 9 *)
{GraphicsRow[{sun2GTransformed, sun2GInverse}], , GraphicsRow[{segmentation[sun2GTransformed], segmentation[sun2GInverse]}]}

(* Input cell 10 *)

 (*Tweaked Segmentation Code using MaxDetect*) 
 
 segmentation12[img_]:=Module[{gaussianFilter, maxValues, centroids, segmented}, gaussianFilter=GaussianFilter[img, 100]; 
 maxValues=MaxDetect[gaussianFilter, Padding -> 1]; 
 centroids=ComponentMeasurements[maxValues, "Centroid"]; 
 segmented=Colorize[WatershedComponents[GradientFilter[img, 2], centroids[[All, 2]]]]] 
 
 GraphicsRow[{segmentation12[sun2GTransformed], segmentation12[sun2GInverse]}] 


(* Input cell 11 *)
(*Tweaked Segmentation Code using MinDetect*) 
 
 segmentation13[img_]:=Module[{gaussianFilter, maxValues, centroids, segmented}, gaussianFilter=GaussianFilter[img, 100]; 
 maxValues=MinDetect[gaussianFilter, Padding -> 1]; 
 centroids=ComponentMeasurements[maxValues, "Centroid"]; 
 segmented=Colorize[WatershedComponents[GradientFilter[img, 2], centroids[[All, 2]]]]] 
 
 GraphicsRow[{segmentation13[sun2GTransformed], segmentation13[sun2GInverse]}]

(* Input cell 12 *)

 
 (*See Section 4.2.1 & Figure 4.3 for Preliminary analysis of Protein Spatial Distribution by Shelley Buchan which was done via vertical and horizontal intensity profiles*) 
 (*The codes below breakdown how the blue and red channels were extracted and how the intensity profiles were drawn*) 
 
 
 (*The allNuclei and allProteins functions below produce a list of the images of all of the individual detected nuclei and their protein distributions, along with their respective centroid coordinates. The input, "blueimg" is the DAPI image, while "redimg" is the Alexa 594 channel image, each from the same protein image dataset. The images are obtained by using the segmentation mask produced by applying the segmentation function detailed above, to the DAPI image in the case of obtaining nuclei images, and the Alexa 594 channel to obtain the protein images.*) 
 
 allNuclei[blueimg_, redimg_]:=allNuclei[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], blueimg}, {"MaskedImage", "Centroid"}]; 
 
 allProteins[blueimg_, redimg_]:=allProteins[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], redimg}, {"MaskedImage", "Centroid"}]; 
 
 (*The imageNuclei and imageProtein functions below ouput an image of only the nucleus/protein selected for analysis (based on the input component number, "number"*) 
 
 imageNucleus[number_, blueimg_, redimg_]:=First@Last@SelectFirst[allNuclei[blueimg, redimg], First[#]==number&];
imageProtein[number_, blueimg_, redimg_]:=First@Last@SelectFirst[allProteins[blueimg, redimg], First[#]==number&]; 
 
 
 (*The componentImageNucleus function below displays an image of the nucleus selected for analysis after using the segmentationOverlay visualisation, while componentImageProtein displays the corresponding protein image, with a white line in the vertical and horizontal axes, respectively. This is done by using the imageProtein and imageNucleus functions defined earlier.*) 
 
 componentImageNucleus[number_, blueimg_, redimg_]:=Show[imageNucleus[number, blueimg, redimg], Graphics[{White, Thick, Line[{{Round[ImageDimensions[imageNucleus[number, blueimg, redimg]][[1]]/2], 0}, {Round[ImageDimensions[imageNucleus[number, blueimg, redimg]][[1]]/2], ImageDimensions[imageNucleus[number, blueimg, redimg]][[2]]}}], Line[{{0, Round[ImageDimensions[imageNucleus[number, blueimg, redimg]][[2]]/2]}, {ImageDimensions[imageNucleus[number, blueimg, redimg]][[1]], Round[ImageDimensions[imageNucleus[number, blueimg, redimg]][[2]]/2]}}]}], ImageSize->Medium]; 
 
 componentImageProtein[number_, blueimg_, redimg_]:=Show[imageProtein[number, blueimg, redimg], Graphics[{White, Thick, Line[{{Round[ImageDimensions[imageProtein[number, blueimg, redimg]][[1]]/2], 0}, {Round[ImageDimensions[imageProtein[number, blueimg, redimg]][[1]]/2], ImageDimensions[imageProtein[number, blueimg, redimg]][[2]]}}], Line[{{0, Round[ImageDimensions[imageProtein[number, blueimg, redimg]][[2]]/2]}, {ImageDimensions[imageProtein[number, blueimg, redimg]][[1]], Round[ImageDimensions[imageProtein[number, blueimg, redimg]][[2]]/2]}}]}], ImageSize->Medium]; 
 
 
 (*The intensitiesComparison functions below plot graphs of the pixel intensities of the selected nucleus (from the intensitiesNucleus function) against the pixel intensities of the corresponding protein distribution (from the intensitiesProtein function) by combining the X functions above for the horizontal (X) profile, and the Y functions above for the vertical (Y) profile.*) 
 
 intensitiesComparisonY[number_, blueimg_, redimg_]:=ListLinePlot[{Transpose[{Range[Length[intensitiesProteinY[number, blueimg, redimg]]], intensitiesProteinY[number, blueimg, redimg]}], Transpose[{Range[Length[intensitiesNucleusY[number, blueimg, redimg]]], intensitiesNucleusY[number, blueimg, redimg]}]}, PlotRange->All, AxesLabel->{"Distance (pixels)", "Grayscale Value"}, PlotStyle->{Red, Blue}, PlotLegends->{"Red Profile", "Blue Profile"}, ImageSize->Medium, AspectRatio->Full, PlotLabel->Style["Vertical Profile of Component "" ToString[number], 20, Bold], PerformanceGoal->"Speed"]; 
 
 intensitiesComparisonX[number_, blueimg_, redimg_]:=ListLinePlot[{Transpose[{Range[Length[intensitiesProteinX[number, blueimg, redimg]]], intensitiesProteinX[number, blueimg, redimg]}], Transpose[{Range[Length[intensitiesNucleusX[number, blueimg, redimg]]], intensitiesNucleusX[number, blueimg, redimg]}]}, PlotRange->All, AxesLabel->{"Distance (pixels)", "Grayscale Value"}, PlotStyle->{Red, Blue}, PlotLegends->{"Red Profile", "Blue Profile"}, ImageSize->Medium, AspectRatio->Full, PlotLabel->Style["Horizontal Profile of Component "" ToString[number], 20, Bold], PerformanceGoal->"Speed"]; 
 
 (*Computing the intensities across the vertical (y-axis) and horizontal (x-axis) nuclei and protein distribution images 
 The intensitiesProtein and intensitiesNucleus functions below produce a list of the pixel intensities of the middle vertical strip, and horizontal strip, of the chosen nucleus/protein image, which have been generated from the functions above. The "Y" in the function name indicates the vertical intensity profile, whlie the "X" indicates the horizontal intensity profile. The images are converted to grayscale to allow comparability of the intensity of the images*) 
 
 intensitiesProteinY[number_, blueimg_, redimg_]:=Transpose[ImageData[ColorConvert[imageProtein[number, blueimg, redimg], "Grayscale"]][[All, All, 1]]][[Round[Length[ImageData[ColorConvert[imageProtein[number, blueimg, redimg], "Grayscale"]][[All, All, 1]]]/2]]]; 
 
 intensitiesNucleusY[number_, blueimg_, redimg_]:=Transpose[ImageData[ColorConvert[imageNucleus[number, blueimg, redimg], "Grayscale"]][[All, All, 1]]][[Round[Length[ImageData[ColorConvert[imageNucleus[number, blueimg, redimg], "Grayscale"]][[All, All, 1]]]/2]]]; 
 
 intensitiesProteinX[number_, blueimg_, redimg_]:=ImageData[ColorConvert[imageProtein[number, blueimg, redimg], "Grayscale"]][[All, All, 1]][[Round[Length[ImageData[ColorConvert[imageProtein[number, blueimg, redimg], "Grayscale"]][[All, All, 1]]]/2]]]; 
 
 intensitiesNucleusX[number_, blueimg_, redimg_]:=ImageData[ColorConvert[imageNucleus[number, blueimg, redimg], "Grayscale"]][[All, All, 1]][[Round[Length[ImageData[ColorConvert[imageNucleus[number, blueimg, redimg], "Grayscale"]][[All, All, 1]]]/2]]]; 
 
 
 (*The plotAndImage function outputs the two graphs generated by the intensitiesComparison functions described earlier, side by side. These plots are followed by the images of the nucleus and protein distributions being analysed (using the componentImage functions described above).
 The plotAndImageMulti function does the same as the plotAndImage function, except it accepts an input of multiple component numbers so that many nuclei of interest can be analysed at once.*) 
 
 plotAndImageMulti[numbers_List, blueimg_, redimg_]:=GraphicsGrid[Table[{intensitiesComparisonX[n, blueimg, redimg], intensitiesComparisonY[n, blueimg, redimg], componentImageNucleus[n, blueimg, redimg], componentImageProtein[n, blueimg, redimg]}, {n, numbers}], Spacings->{1, 1}]; 
 
 
 plotAndImageMulti[{149, 70}, sun2Blue, sun2Red] (*Plotting intensity profiles of extrafusal and intrafusal nuclei and proteins corresponding to Figure 4.3*) 
 


(* Input cell 13 *)
(*In Section 4, The MNIST Digit Dataset was used as a benchmark to validate computational pipeline (Figure 4.1) 
 loaded and filtered to only obtain 0, 1 and 6. see below codes *) 
 
 
 imgs= ResourceData["MNIST"]; (*loads MNIST Dataset*) 
 Dimensions[imgs] (*checks Dimensions*) 
 
 (*Filtering to Obtain only 0, 1 and 6 images with the function "Keys"*) 
 zeroes=Select[imgs, # [[2]]==0&]; 
 ones = Select[imgs, # [[2]]==1&]; 
 sixes=Select[imgs, # [[2]]==6&]; 
 zeroesK =Keys[zeroes]; 
 onesK=Keys[ones]; 
 sixesK = Keys[sixes]; 
 
 (* Loading & Truncating a simple neural network, LeNet (see Section 4.4.1) to act as feature extraction network instead of a classifier*) 
 
 lenet=NetModel["LeNet"]; 
 lenetIn = NetInitialize[lenet]; (*adds random weights and biases*) 
 lenetF = Take[lenetIn, {1, -5}]; (*takes all layers from the top to the "flatten layer", which is the feature extraction layer*) 
 
 (* Apply truncated LeNet to all images above (those obtained with the function "Keys") to obtain feature vectors and merge all vectors in pairs (ie. "0 vs 1" and "0 vs 6" as described in Section 4.3 and Sections 5.1.1.1 and 5.1.1.2 In code below, Z, 1 and 6 which end features correspond to zero, one and 6 digits respectively *) 
 
 featuresZ=lenetF/@zeroesK; 
 features1=lenetF/@onesK; 
 features6 = lenetF/@sixesK; 
 Length[featuresZ] (*check the Length of features zero (0)*) 
 Length[features1] (*check the Length of features one (1)*) 
 Length[features6] (*check the Length of features six (6)*) 


(* Input cell 14 *)
(*Merging the features 0 vs 1 and 0 vs 6 ; and checking dimensions too*) 
 allfeaturesZ1=Join[featuresZ, features1]; 
 allfeaturesZ6 = Join[featuresZ, features6]; 
 Dimensions[allfeaturesZ1] (*check the Dimensions of feature digit pairs 0 vs 1*) 
 Dimensions[allfeaturesZ6] (*check the Dimensions of feature digit pairs 0 vs 6*)

(* Input cell 15 *)

 (*Codes below correspond to Dimensionality reduction in Section 4.5 *) 
 
 redMZ1=DimensionReduce[allfeaturesZ1, 2, Method->"TSNE"]; (*reduction in dimensions for pair 0 vs 1*) 
 redMZ6=DimensionReduce[allfeaturesZ6, 2, Method->"TSNE"]; (*reduction in dimensions for pair 0 vs 6*) 
 
 (*Checking dimensions of the reduced paired features*) 
 
 Dimensions[redMZ1] 
 Dimensions[redMZ6] 


(* Input cell 16 *)
(*Defining a funtion block that does:
 
 a) Automatic Clustering (based on k you specify, but k is set to 2 initially), b) Plotting of reduced tSNE features, c) Calculation of Confusion Matrix for MNIST Digit Pairs 
 *) 
 
 Clear[ClusterAndPlotWithConfusion]

 ClusterAndPlotWithConfusion[featuresA_, featuresB_, k_:2, opts:OptionsPattern[{"Method"-"KMeans"}], labels_List:{"Image 0", "Image 1"}]:=Module[{allFeatures, reduced, clusters, clusA, clusB, cmatrix, rowLabels, colLabels, rowTotals, maxIndices, colorMap, confusionGrid, scatterPlot, finalLayout, statsText, headerLabels}, (*Step 1:Join feature vectors*) 
 allFeatures=Join[featuresA, featuresB]; 
 
 (*Step 2:Dimensionality reduction with t-SNE*) 
 reduced=DimensionReduce[allFeatures, 2, Method->"TSNE"]; 
 
 (*Step 3:Clustering*) 
 clusters=ClusteringComponents[allFeatures, k, 1, opts]; 
 
 (*Step 4:Separate clusters*) 
 clusA=Take[clusters, Length[featuresA]]; 
 clusB=Take[clusters, -Length[featuresB]]; 
 
 (*Step 5:Confusion matrix*) 
 cmatrix={Table[Count[clusA, j], {j, 1, k}], Table[Count[clusB, j], {j, 1, k}]}; 
 rowLabels=labels; 
 rowTotals=Total[cmatrix, {2}]; 
 maxIndices=Ordering[#, -1]&/@cmatrix; 
 colorMap=Association[Table[i->ColorData[97][i], {i, k}]]; 
 
 (*Step 6:Confusion matrix column headers*) 
 colLabels=Table[Grid[{{"Predicted Cluster ""ToString[i]}}, Alignment->Center], {i, 1, k}]; 
 headerLabels=Join[{""}, colLabels, {Grid[{{"Total Samples Per Image"}}, Alignment->Center], "Accuracy (%)", "Error (%)"}]; 
 
 (*Step 7:Confusion matrix grid with shading and formatting*) 
 confusionGrid=Grid[Prepend[Table[With[{row=cmatrix[[i]], total=rowTotals[[i]], maxIndex=First@maxIndices[[i]], label=rowLabels[[i]]}, Module[{rowStyled, acc, err}, acc=N[100*row[[maxIndex]]/total, 3]; 
 err=N[100*(1-row[[maxIndex]]/total), 3]; 
 rowStyled=Table[If[j==maxIndex, Item[Style[row[[j]], Bold, FontFamily->"Arial", FontSize->16], Background->Lighter[colorMap[j], 0.9], Alignment->Center, ItemSize->{8, 2}], Item[Style[row[[j]], FontFamily->"Arial", FontSize->16], Alignment->Center, ItemSize->{8, 2}]], {j, Length[row]}]; 
 Join[{Style[label, Bold, FontFamily->"Arial", FontSize->16]}, rowStyled, {Style[total, Italic, Gray, FontSize->15], Style[acc, Bold, Italic, Black, FontSize->15], Style[err, Bold, Italic, Black, FontSize->15]}]]], {i, Length[cmatrix]}], headerLabels], Frame->All, Alignment->Center, Spacings->{2, 1.5}, ItemSize->{{10, 10, 10, 10, 10, 10}, 2.5}, BaseStyle->{FontFamily->"Arial", FontSize->16}]; 
 
 (*Step 8:Unified stats+legend box*) statsText=Column[{Row[{Style["Length: ", Bold], ToString[Length[clusters]]}], Row[{Style["Dimensions: ", Bold], ToString[Dimensions[clusters]]}], Spacer[10], Style["Legend:", Bold], Column[Table[Row[{Style["[FilledCircle]", FontColor->colorMap[i], FontSize->12], " Cluster ""ToString[i]}], {i, 1, k}]]}, Alignment->Left, Spacings->0.5]; 
 
 (*Step 9: Scatter plot with aligned statsText*) scatterPlot=Show[ListPlot[MapThread[Style[#1, colorMap[#2]]&, {reduced, clusters}], PlotStyle->PointSize[Medium], AxesLabel->{"tSNE-1", "tSNE-2"}, ImageSize->500, PlotRangeClipping->False, ImagePadding->{{50, 200}, {50, 50}}], Epilog->{Inset[Style[statsText, FontSize->12, FontFamily->"Arial"], {60, -55}, {Left, Bottom} (*alignment*)]}]; 
 
 
 (*Step 10: Combine plot and confusion matrix with divider*) finalLayout=Grid[{{scatterPlot, Column[{Style["Confusion Matrix", Bold, 18], confusionGrid}, Alignment->Center]}}, Alignment->Center, Spacings->{5, 2}, Dividers->{{False, Directive[Thick, Dashed], False}, {False}}]; 
 
 
 (*Step 11:Return plot and confusion matrix separately*) 
 
 Return[{scatterPlot, confusionGrid}]]

(* Input cell 17 *)

 (*Clustering with Euclidean Distance, Displaying tSNE reduced scatter plot and plotting confusion matrix for Digit Pairs. Corresponds to Section 5.1.1.1 and 5.1.1.2. The Figures can also be seen in Figures 5.1 and 5.3 as well as Tables 5.1 and 5.3*) 
 
 {plot2a, table2a}=ClusterAndPlotWithConfusion[featuresZ, features1, 2, "Method"-"KMeans", DistanceFunction->EuclideanDistance, {"Image 0", "Image 1"}]; 
 {plot2b, table2b}=ClusterAndPlotWithConfusion[featuresZ, features6, 2, "Method"-"KMeans", DistanceFunction->EuclideanDistance, {"Image 0", "Image 6"}]; 
 
 (*Defining a function that displays scatter plot with a title you define*) 
 addSubtitle[plot_, subtitle_]:=
 Column[{Style[subtitle, Bold, 16, Black, FontFamily->"Arial"], plot}, Alignment->Center, Spacings->1]; 
 
 addSubtitle[plot2a, "<Clustering Digits 0 & 1 into 2 (KMeans, Level=1, Euclidean 
Distance)"] 
 addSubtitle[plot2b, "<Clustering Digits 0 & 6 into 2 (KMeans, Level=1, Euclidean 
Distance)"] 
 
 (*Displaying Confusion tables as well*) 
 table2a
table2b 


(* Input cell 18 *)
(*Section corresponds to Supplementary Material Figures S2 and S4*) 
 (* Clustering with Default Mathematica Distance Metric, Displaying tSNE reduced scatter plot and plotting confusion matrix for Digit Pairs *) 
 
 {plot1a, table1a}=ClusterAndPlotWithConfusion[featuresZ, features1, 2, "Method"-"KMeans", {"Image 0", "Image 1"}]; 
 {plot1b, table1b}=ClusterAndPlotWithConfusion[featuresZ, features6, 2, "Method"-"KMeans", {"Image 0", "Image 6"}]; 
 
 addSubtitle[plot1a, "Clustering Digits 0 & 1 into 2 (KMeans, Level=1)"] 
 addSubtitle[plot1b, "Clustering Digits 0 & 6 into 2 (KMeans, Level=1)"] 
 
 table1a 
 table1b

(* Input cell 19 *)
(*Using the loaded MNIST Digit Dataset above in Section 4 above, a more sophisticated network, ResNet-50 (See Section 4.4.2), was used for the same pipeline above (See Figure 4.1). See below codes *) 
 
 imgs= ResourceData["MNIST"]; 
 zeroes=Select[imgs, # [[2]]==0&]; 
 ones = Select[imgs, # [[2]]==1&]; 
 sixes=Select[imgs, # [[2]]==6&]; 
 zeroesK =Keys[zeroes]; 
 onesK=Keys[ones]; 
 sixesK = Keys[sixes];

(* Input cell 20 *)

 (* Loading & Truncating a complex neural network, ResNet (see Section 4.4.2) to act as feature extraction network instead of a classifier*) 
 
 resnet=NetModel["ResNet-50 Trained on ImageNet Competition Data"]; 
 resnetF=Take[resnet, {1, -3}]; (*takes all layers from the top to the "flatten layer", which is the feature extraction layer*) 
 
 (* Apply truncated ResNet to all images above (those obtained with the function "Keys") to obtain feature vectors and merge all vectors in pairs (ie. "0 vs 1" and "0 vs 6" as described in Section 4.3 and Sections 5.1.1.1 and 5.1.1.2 In code below, Z, 1 and 6 which end features correspond to zero, one and 6 digits respectively *) 
 
 rfeatZ=resnetF/@zeroesK; 
 rfeat1=resnetF/@onesK; 
 rfeat6 = resnetF/@sixesK; 


(* Input cell 21 *)

 (*For ResNet-50*) 
 (*Using predefined functions "ClusterAndPlotWithConfusion" and addSubtitle below*) 
 (*Clustering with Euclidean Distance, Displaying tSNE reduced scatter plot and plotting confusion matrix for Digit Pairs. Corresponds to Section 5.1.1.1 and 5.1.1.2. The Figures can also be seen in Figures 5.2 and 5.4 as well as Tables 5.2 and 5.4*) 
 
 {plot2Ra, table2Ra}=ClusterAndPlotWithConfusion[rfeatZ, rfeat1, 2, "Method"-"KMeans", DistanceFunction->EuclideanDistance, {"Image 0", "Image 1"}]; 
 {plot2Rb, table2Rb}=ClusterAndPlotWithConfusion[rfeatZ, rfeat6, 2, "Method"-"KMeans", DistanceFunction->EuclideanDistance, {"Image 0", "Image 6"}]; 
 
 addSubtitle[plot2Ra, "<Clustering Digits 0 & 1 into 2 (KMeans, Level=1, Euclidean 
Distance)"]
addSubtitle[plot2Rb, "<Clustering Digits 0 & 6 into 2 (KMeans, Level=1, Euclidean 
Distance)"] 
 
 (*Displaying Confusion tables as well*) 
 table2Ra 
 table2Rb 


(* Input cell 22 *)

 (*Section corresponds to Supplementary Material Figures S3 and S5*) 
 (* Clustering with Default Mathematica Distance Metric, Displaying tSNE reduced scatter plot and plotting confusion matrix for Digit Pairs *) 
 
 {plot1Ra, table1Ra}=ClusterAndPlotWithConfusion[rfeatZ, rfeat1, 2, "Method"-"KMeans", {"Image 0", "Image 1"}]; 
 {plot1Rb, table1Rb}=ClusterAndPlotWithConfusion[rfeatZ, rfeat6, 2, "Method"-"KMeans", {"Image 0", "Image 6"}]; 
 
 addSubtitle[plot1Ra, "Clustering Digits 0 & 1 into 2 (KMeans, Level=1)"]
addSubtitle[plot1Rb, "Clustering Digits 0 & 6 into 2 (KMeans, Level=1)"] 
 
 (*Displaying Confusion tables as well*) 
 table2Ra 
 table2Rb 


(* Input cell 23 *)

 resnet=NetModel["ResNet-50 Trained on ImageNet Competition Data"]; 
 resnetF=Take[resnet, {1, -3}]; 
 
 (*Extracting Proteins*) 
 allProteins[blueimg_, redimg_]:=allProteins[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], redimg}, {"MaskedImage", "Centroid"}];
prolist=allProteins[sun2Blue, sun2Red]; 
 prolistV=Values[prolist]; (*Takes only the values containing images and centroid co-ordinates*) 
 proPics=prolistV[[All, 1]]; (*Pick only protein images*) 
 
 rProtein=resnetF/@proPics; 
 Dimensions[rProtein] 


(* Input cell 24 *)

 
 (*Defining a funtion block called "ProteinClusterAndPlot" for Proteins that does:
 
 a) Automatic Clustering (based on k you specify, but k is set to 2 initially), b) Plotting of reduced tSNE features, c) Returns only the scatter plots and clusters for proteins only 
 *)

(* Input cell 25 *)
ProteinClusterAndPlot[proteinFeatures_, k_:2, opts:OptionsPattern[]]:=Module[{reduced, clusters, tally, statsText, plot}, (*Step 1:Reduce features with t-SNE*) 
 reduced=DimensionReduce[proteinFeatures, 2, Method->"TSNE"]; 
 
 (*Step 2:Clustering*) 
 clusters=ClusteringComponents[proteinFeatures, k, 1, opts]; 
 
 (*Step 3:Count proteins per cluster*) 
 tally=Tally[clusters]; 
 
 (*Step 4:Stats block*) statsText=Column[{Row[{Style["Length: ", Bold], ToString[Length[clusters]]}], Row[{Style["Dimensions: ", Bold], ToString[Dimensions[clusters]]}], "", Style["Protein Clusters:", Bold], Column[Map[Row[{" Cluster ", # [[1]], " [RightArrow] ", Style[# [[2]], Bold], " proteins"}]&, tally]], "", Style["Legend:", Bold], Column[Table[Row[{" ", Style["[FilledCircle]", FontColor->ColorData[97][i], FontSize->14], " Cluster ""ToString[i]}], {i, 1, k}]]}, Alignment->Left, Spacings->0.6]; 
 
 (*Step 5:Scatter plot with stats positioned further to the right*) plot=ListPlot[MapThread[Style[#1, ColorData[97][#2]]&, {reduced, clusters}], PlotStyle->PointSize[Medium], AxesLabel->{"tSNE-1", "tSNE-2"}, PlotRange->Automatic, ImageSize->500, ImagePadding->{{50, 180}, {35, 40}}, (*extra room on right for stats*) PlotRangeClipping->False, Epilog->{Inset[Style[statsText, FontSize->12, FontFamily->"Arial"], Scaled[{1.05, 0.05}], (*shift further right& slightly upward*) {Left, Bottom}]}]; 
 
 (*Step 6:Return both plot and clusters*) 
 {plot, clusters}]

(* Input cell 26 *)

 (*Clustering Proteins into 2 or 3 cluster with Euclidean Distance. See Section 5 and Figures 5, 6*) 
 
 {myPlot1a, clus1a}=ProteinClusterAndPlot[rProtein, 2, Method->"KMeans", DistanceFunction->EuclideanDistance]; 
 {myPlot1b, clus1b}=ProteinClusterAndPlot[rProtein, 3, Method->"KMeans", DistanceFunction->EuclideanDistance]; 
 
 (*Displaying scatter plots*) 
 Column[{Style["<ResNet-50 Feature Vector Visualizations of SUN2 Proteins 
Only In Clusters", Black, Bold, 20, FontFamily->"Arial"], addSubtitle[myPlot1a, "<Clustering into 2 (KMeans, Level=1, Euclidean 
Distance)"], addSubtitle[myPlot1b, "<Clustering into 3 (KMeans, Level=1, Euclidean 
Distance)"]}, Spacings->2, Alignment->Center] 


(* Input cell 27 *)
(*Overlaying Clusters coded with Colours on actual SUN2Brightfield Images*) 
 
 segmentationOverlay[blueimg_, redimg_, greyimg_]:=Module[{segment=segmentation[blueimg], measurements}, measurements=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segment], 10], 5], redimg}, {"Centroid", "MaskedImage"}]; (*The variable "measurements" is for the tooltips function, so that the tooltip will display the component number and the corresponding protein image when you hover over a centroid*) 
 HighlightImage[greyimg+MorphologicalPerimeter[segment], Graphics[{Blue, PointSize[0.005], Table[Tooltip[Style[Point[measurements[[i, 2, 1]]], PointSize[0.008]], Column[{Style["Component ""ToString[i], Bold], measurements[[i, 2, 2]]}]], {i, Length[measurements]}]}]]]; (*Function here overlays centroids without color (ie. uses dots) and gives tooltip showing segment and its component number*) 


(* Input cell 28 *)
(*Defined function from the overlay function above that works with any number of clusters (k) specified and shows color coded clusters on overlays*) 
 
 segmentationOverlayWithClustersUniversal[blueimg_, redimg_, greyimg_, clusterLabels_]:=Module[{segment, measurements, numClusters, colors}, (*Step 1:Segmentation and component measurements*) segment=segmentation[blueimg]; 
 measurements=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segment], 10], 5], redimg}, {"Centroid", "MaskedImage"}]; 
 (*Step 2:Detect number of unique clusters*) numClusters=Max[clusterLabels]; 
 (*Step 3:Use the same shifted ColorData[97] as t-SNE (#2+1 rule)*) colors=Table[ColorData[97][k], {k, 1, numClusters}]; 
 (*Step 4:Overlay with matching colors+tooltip*) HighlightImage[greyimg+MorphologicalPerimeter[segment], Graphics[Table[{colors[[clusterLabels[[i]]]], PointSize[0.008], Tooltip[Style[Point[measurements[[i, 2, 1]]], PointSize[0.008]], Column[{Style["Component ""ToString[i], Bold], measurements[[i, 2, 2]]}]]}, {i, Length[measurements]}]]]] 


(* Input cell 29 *)
(*Displaying overlays of SUN2 Protein only clusters on TRANS image *) 
 
 Grid[{{Style["<Overlay with 2 Clusters of Proteins Only (KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Proteins Only (KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[sun2Blue, sun2Red, sun2Grey, clus1a], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[sun2Blue, sun2Red, sun2Grey, clus1b], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 30 *)
(*Extracting Protein Segments*) 
 allProteins[blueimg_, redimg_]:=allProteins[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], redimg}, {"MaskedImage", "Centroid"}];
prolist=allProteins[sun2Blue, sun2Red]; 
 prolistV=Values[prolist]; 
 proPics=prolistV[[All, 1]]; 
 
 (*Extracting Nuclei Segments*) 
 allNuclei[blueimg_, redimg_]:=allNuclei[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], blueimg}, {"MaskedImage", "Centroid"}]; 
 nucList=allNuclei[sun2Blue, sun2Red]; 
 nuPics = Values[nucList][[All, 1]]; 
 
 (*Pairing them side by side for each component to obtain 4096 feature vectors instead of 2048 for protein only vectors.*) 
 allNuPro=Transpose[{nuPics, proPics}]; 
 rNuPro=resnetF/@allNuPro; 
 Dimensions[rNuPro] 
 
 mergedFeatures=Map[Join[# [[1]], # [[2]]]&, (*take nucleus (first) then protein (second)*) 
 rNuPro]; 
 
 Dimensions[mergedFeatures]

(* Input cell 31 *)

 (*Defining a function for Combined Segments that does same thing as "ProteinClusterAndPlot"*) 
 
 NucleusProteinClusterAndPlot[nucleusProteinFeatures_, k_:2, opts:OptionsPattern[]]:=Module[{reduced, clusters, tally, statsText, plot}, (*Step 1:Reduce features with t-SNE*) 
 reduced=DimensionReduce[nucleusProteinFeatures, 2, Method->"TSNE"]; 
 
 (*Step 2:Clustering*) 
 clusters=ClusteringComponents[nucleusProteinFeatures, k, 1, opts]; 
 
 (*Step 3:Count proteins per cluster*) 
 tally=Tally[clusters]; 
 
 (*Step 4:Stats block*) statsText=Column[{Row[{Style["Length: ", Bold], ToString[Length[clusters]]}], Row[{Style["Dimensions: ", Bold], ToString[Dimensions[clusters]]}], "", Style["Nuclei+Protein Clusters:", Bold], Column[Map[Row[{" Cluster ", # [[1]], " [RightArrow] ", Style[# [[2]], Bold], " proteins"}]&, tally]], "", Style["Legend:", Bold], Column[Table[Row[{" ", Style["[FilledCircle]", FontColor->ColorData[97][i], FontSize->14], " Cluster ""ToString[i]}], {i, 1, k}]]}, Alignment->Left, Spacings->0.6]; 
 
 (*Step 5:Scatter plot with stats positioned further to the right*) plot=ListPlot[MapThread[Style[#1, ColorData[97][#2]]&, {reduced, clusters}], PlotStyle->PointSize[Medium], AxesLabel->{"tSNE-1", "tSNE-2"}, PlotRange->Automatic, ImageSize->500, ImagePadding->{{50, 180}, {35, 40}}, (*extra room on right for stats*) PlotRangeClipping->False, Epilog->{Inset[Style[statsText, FontSize->12, FontFamily->"Arial"], Scaled[{1.05, 0.05}], (*shift further right& slightly upward*) {Left, Bottom}]}]; 
 
 (*Step 6:Return both plot and clusters*) 
 {plot, clusters}]

(* Input cell 32 *)
{{myPlot2a, clus2a}=NucleusProteinClusterAndPlot[mergedFeatures, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {myPlot2b, clus2b}=NucleusProteinClusterAndPlot[mergedFeatures, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];, , Column[{Style["<ResNet-50 Feature Vector Visualizations of SUN2 Nuclei + 
Proteins In Clusters", Black, Bold, 20, FontFamily->"Arial"], addSubtitle[myPlot2a, "<Clustering into 2 (KMeans, Level=1, Euclidean 
Distance)"], addSubtitle[myPlot2b, "<Clustering into 3 (KMeans, Level=1, Euclidean 
Distance)"]}, Spacings->2, Alignment->Center]}

(* Input cell 33 *)
Grid[{{Style["<Overlay with 2 Clusters of Nuclei & Proteins (KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Nuclei & Proteins (KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[sun2Blue, sun2Red, sun2Grey, clus2a], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[sun2Blue, sun2Red, sun2Grey, clus2b], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 34 *)
(*Considering the fact that all segments have different dimensions and they were in the RGBA color space, there was the need to streamline all inputs to ResNet-50 to meet its input specifications for better results*)

(* Input cell 35 *)
{ResizeImageForResNet[img_]:=ColorConvert[ImageResize[RemoveAlphaChannel[img], {224, 224}], "RGB"];, , nucleusResized=ResizeImageForResNet/@nuPics;, , proteinResized=ResizeImageForResNet/@proPics;, , Dimensions[ImageData[nucleusResized[[5]]]]}

(* Input cell 36 *)
{ExtractFeatures[img_]:=resnetF[img];, , nucleusFeatures=ExtractFeatures/@nucleusResized;, , proteinFeatures=ExtractFeatures/@proteinResized;, }

(* Input cell 37 *)
{featureVectors=MapThread[Join, {nucleusFeatures, proteinFeatures}];, , Dimensions[featureVectors], }

(* Input cell 38 *)
{{myNPlot, NproteinClusters}=NucleusProteinClusterAndPlot[featureVectors, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {myNPlot2, NproteinClusters2}=NucleusProteinClusterAndPlot[featureVectors, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];}

(* Input cell 39 *)
Column[{Style["<ResNet-50 Feature Vector Visualizations of Resized SUN2 
Nuclei + Proteins In Clusters", Black, Bold, 20, FontFamily->"Arial"], addSubtitle[myNPlot, "<Clustering into 2 (KMeans, Level=1, Euclidean 
Distance)"], addSubtitle[myNPlot2, "<Clustering into 3 (KMeans, Level=1, Euclidean 
Distance)"]}, Spacings->2, Alignment->Center]

(* Input cell 40 *)
Grid[{{Style["<Overlay with 2 Clusters of Resized Nuclei & Proteins 
(KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Resized Nuclei & Proteins 
(KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[sun2Blue, sun2Red, sun2Grey, NproteinClusters], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[sun2Blue, sun2Red, sun2Grey, NproteinClusters2], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 41 *)
(* =========================1) 
 Correct resize for ResNet-50 (aspect ratio preserved+padded)=========================*) 
 ClearAll[ResizeCorrect];

 ResizeCorrect[img_, size_:224]:=Module[{im, dims, scale, newDims, resized, pad1, pad2}, im=RemoveAlphaChannel@ColorConvert[img, "RGB"]; 
 dims=ImageDimensions[im]; 
 scale=If[Max[dims]>size, N[size/Max[dims]], 1.0]; 
 newDims=Round[scale*dims]; 
 newDims=Map[Max[#, 1]&, newDims]; 
 resized=If[scale==1.0, im, ImageResize[im, newDims, Resampling->"Lanczos"]]; 
 pad1=Floor[(size-newDims)/2]; 
 pad2=(size-newDims)-pad1; 
 ImagePad[resized, {{pad1[[1]], pad2[[1]]}, {pad1[[2]], pad2[[2]]}}, 0]]; 
 
 (* =========================2) 
 WRONG method (forced squish to 224 x224)=========================*) 
 
 ClearAll[ResizeWrong];

 ResizeWrong[img_]:=ColorConvert[ImageResize[RemoveAlphaChannel[img], {224, 224}], "RGB"]; 
 
 (* =========================3) 
 Comparison table with "Wrong" and "Correct" in the HEADLINE=========================*) 

ClearAll[CompareResizeTable];

 CompareResizeTable[img_]:=Module[{original, wrong, correct}, original=img; 
 wrong=ResizeWrong[original]; 
 correct=ResizeCorrect[original]; 
 Grid[{{"", Style["Original", Bold, 14], Style["WRONG: Forced resize to 224[Times]224", Red, Bold, 14], Style["CORRECT: Aspect ratio preserved + padded", Darker[Green], Bold, 14]}, {"Dimensions", Dimensions[ImageData[original]], Dimensions[ImageData[wrong]], Dimensions[ImageData[correct]]}, {"Image", original, wrong, correct}}, Frame->All, Alignment->Center, ItemStyle->Directive[12]]];

(* Input cell 42 *)
{i=5;, , CompareResizeTable[nuPics[[i]]], }

(* Input cell 43 *)
{nucleusResized=ResizeCorrect/@nuPics;, , proteinResized=ResizeCorrect/@proPics;, }

(* Input cell 44 *)
nuPics[[;; 5]]

(* Input cell 45 *)
nucleusResized[[;; 5]]

(* Input cell 46 *)
proteinResized[[;; 5]]

(* Input cell 47 *)
Manipulate[Module[{original, wrong, correct}, original=nuPics[[i]]; 
 wrong=ResizeWrong[original]; 
 correct=ResizeCorrect[original]; 
 Grid[{{"", Style["Original", Bold, 14], Style["WRONG: Forced 224[Times]224", Red, Bold, 14], Style["CORRECT: Preserve + Pad", Darker[Green], Bold, 14]}, {"Index", i, "", ""}, {"Dimensions", Dimensions[ImageData[original]], Dimensions[ImageData[wrong]], Dimensions[ImageData[correct]]}, {"Image", original, wrong, correct}}, Frame->All, Alignment->Center, ItemStyle->Directive[12], Spacings->{1.2, 1.2}]], {{i, 1, "Image #"}, 1, Length[nuPics], 1, Appearance->"Labeled"}]

(* Input cell 48 *)
{ClearAll[BlockRows, ContinuousManipulateStyleTable];, , BlockRows[i_, thumb_:140]:=Module[{orig, wrong, corr}, orig=nuPics[[i]]; 
 wrong=ResizeWrong[orig]; 
 corr=ResizeCorrect[orig]; 
 {{"Index", i, "", ""}, {"Dimensions", Dimensions[ImageData[orig]], Dimensions[ImageData[wrong]], Dimensions[ImageData[corr]]}, {"Image", Show[orig, ImageSize->thumb], Show[wrong, ImageSize->thumb], Show[corr, ImageSize->thumb]}}];, , ContinuousManipulateStyleTable[start_:1, nImages_:5, thumb_:140]:=Module[{stop, idx, rows}, stop=Min[start+nImages-1, Length[nuPics]]; 
 idx=Range[start, stop]; 
 rows=Join@@(BlockRows[#, thumb]&/@idx); 
 Grid[Join[{{"", Style["Original", Bold, 14], Style["WRONG: Forced 224[Times]224", Red, Bold, 14], Style["CORRECT: Preserve + Pad", Darker[Green], Bold, 14]}}, rows], Frame->All, Dividers->All, (*gives you border lines everywhere*) Alignment->Center, Spacings->{0.8, 0.6}, (*tighten horizontal/vertical spacing*) ItemSize->{{8, 18, 18, 18}, Automatic}, (*keep columns consistent*) FrameStyle->Directive[GrayLevel[0.6], Thickness[0.002]], BaseStyle->{12}]];, , Manipulate[Module[{start}, start=(page-1)*nPerPage+1; 
 ContinuousManipulateStyleTable[start, nPerPage, thumb]], {{page, 1, "Page"}, 1, Dynamic[Ceiling[Length[nuPics]/nPerPage]], 1, Appearance->"Labeled"}, {{nPerPage, 5, "Images per page"}, {5, 10, 15}, ControlType->PopupMenu}, {{thumb, 140, "Thumb size"}, {110, 140, 170}, ControlType->PopupMenu}], }

(* Input cell 49 *)
(*Adds clear border lines for each "lane" (row block)+full cell borders.Uses your existing:ResizeWrong, ResizeCorrect, nuPics*) ClearAll[BlockRows, ContinuousManipulateStyleTable];

 BlockRows[i_, thumb_:140]:=Module[{orig, wrong, corr}, orig=nuPics[[i]]; 
 wrong=ResizeWrong[orig]; 
 corr=ResizeCorrect[orig]; 
 {{"Index", i, "", ""}, {"Dimensions", Dimensions[ImageData[orig]], Dimensions[ImageData[wrong]], Dimensions[ImageData[corr]]}, {"Image", Show[orig, ImageSize->thumb], Show[wrong, ImageSize->thumb], Show[corr, ImageSize->thumb]}}];

 ContinuousManipulateStyleTable[start_:1, nImages_:5, thumb_:140]:=Module[{stop, idx, rows, laneBreakRows}, stop=Min[start+nImages-1, Length[nuPics]]; 
 idx=Range[start, stop]; 
 rows=Join@@(BlockRows[#, thumb]&/@idx); 
 (*lane breaks:header row is 1, then every 3 rows after that ends a lane*) laneBreakRows=Join[{1}, 1+3 Range[nImages]]; 
 Grid[Join[{{"", Style["Original", Bold, 14], Style["WRONG: Forced 224[Times]224", Red, Bold, 14], Style["CORRECT: Preserve + Pad", Darker[Green], Bold, 14]}}, rows], Frame->All, Dividers->{All, (*vertical lines for all columns*) Join[{All}, (*default thin lines between rows*) Thread[laneBreakRows->Directive[GrayLevel[0.2], Thickness[0.004]]] (*thick lane lines*)]}, Alignment->Center, Spacings->{0.8, 0.6}, ItemSize->{{10, 18, 18, 18}, Automatic}, FrameStyle->Directive[GrayLevel[0.5], Thickness[0.003]], BaseStyle->{12}]];

 Manipulate[Module[{start}, start=(page-1)*nPerPage+1; 
 ContinuousManipulateStyleTable[start, nPerPage, thumb]], {{page, 1, "Page"}, 1, Dynamic[Ceiling[Length[nuPics]/nPerPage]], 1, Appearance->"Labeled"}, {{nPerPage, 5, "Images per page"}, {5, 10, 15}, ControlType->PopupMenu}, {{thumb, 140, "Thumb size"}, {110, 140, 170}, ControlType->PopupMenu}]

(* Input cell 50 *)
(* ================================COMPONENT[RightArrow] CLUSTER TABLE Proteins only, Proteins+Nuclei, Truncated ResNet================================*) 
 (*Cluster label vectors already exist:clus1a, clus1b, clus2a, clus2b, NproteinClusters, NproteinClusters2*) labels={clus1a, (*ResNet-50 Proteins only, k=2*) clus1b, (*ResNet-50 Proteins only, k=3*) 
 clus2a, (*ResNet-50 Nuclei+Proteins, k=2*) 
 clus2b, (*ResNet-50 Nuclei+Proteins, k=3*) 
 NproteinClusters, (*Truncated ResNet-50, k=2*) 
 NproteinClusters2 
 (*Truncated ResNet-50, k=3*)};

 headers={"Component ID", "Proteins (ResNet-50, k = 2)", "Proteins (ResNet-50, k = 3)", "Nuclei + Proteins (ResNet-50, k = 2)", "Nuclei + Proteins (ResNet-50, k = 3)", "Proteins (Truncated ResNet-50, k = 2)", "Proteins (Truncated ResNet-50, k = 3)"};

 table=Prepend[Table[Prepend[labels[[All, i]], i], {i, Length[labels[[1]]]}], headers];

 Grid[table, Frame->All, Alignment->Center, ItemStyle->Directive[FontFamily->"Arial", 11], Spacings->{2, 1.2}]

(* Input cell 51 *)
(* ================================COMPONENT[RightArrow] CLUSTER TABLE (+Dominant Cluster) Proteins only, Proteins+Nuclei, Truncated ResNet================================*) labels={clus1a, (*ResNet-50 Proteins only, k=2*) clus1b, (*ResNet-50 Proteins only, k=3*) clus2a, (*ResNet-50 Nuclei+Proteins, k=2*) clus2b, (*ResNet-50 Nuclei+Proteins, k=3*) NproteinClusters, (*Truncated ResNet-50, k=2*) NproteinClusters2 (*Truncated ResNet-50, k=3*)};

 headers={"Component ID", "Proteins (ResNet-50, k = 2)", "Proteins (ResNet-50, k = 3)", "Nuclei + Proteins (ResNet-50, k = 2)", "Nuclei + Proteins (ResNet-50, k = 3)", "Proteins (Truncated ResNet-50, k = 2)", "Proteins (Truncated ResNet-50, k = 3)", "Dominant Cluster"}; 

(*Dominant cluster for one component=most frequent label across methods*) 
 dominantCluster[row_List]:=Module[{t=Tally[row], max}, max=Max[t[[All, 2]]]; 
 (*If tie:choose the smallest cluster number among the tied ones*) Min[Cases[t, {c_, n_}/; n==max:>c]]];

 table=Prepend[Table[Module[{row=labels[[All, i]]}, Join[{i}, row, {dominantCluster[row]}]], {i, Length[labels[[1]]]}], headers];

 Grid[table, Frame->All, Alignment->Center, ItemStyle->Directive[FontFamily->"Arial", 11], Spacings->{2, 1.2}]

(* Input cell 52 *)
(* ================================COMPONENT[RightArrow] CLUSTER TABLE (+Dominant Cluster)================================*) labels={clus1a, clus1b, clus2a, clus2b, NproteinClusters, NproteinClusters2};

 headers={"Component ID", "Proteins (ResNet-50, k = 2)", "Proteins (ResNet-50, k = 3)", "Nuclei + Proteins (ResNet-50, k = 2)", "Nuclei + Proteins (ResNet-50, k = 3)", "Proteins (Truncated ResNet-50, k = 2)", "Proteins (Truncated ResNet-50, k = 3)", "Dominant Cluster"}; 

(*Dominant cluster=most frequent label per component*) 
 dominantCluster[row_List]:=Module[{t=Tally[row], m}, m=Max[t[[All, 2]]]; 
 Min[Cases[t, {c_, n_}/; n==m:>c]]];

 table=Prepend[Table[Module[{row=labels[[All, i]]}, Join[{i}, row, {dominantCluster[row]}]], {i, Length[labels[[1]]]}], headers];

 Grid[table, Frame->All, Alignment->Center, ItemStyle->{{1, All}->Bold, (*header row bold*) {All, -1}->Bold (*Dominant Cluster column bold*)}]

(* Input cell 53 *)
(*Replicating steps for spatial distribution above on 
 a) Lamin-B1 proteins only 
 b) Combined Lamin-B1 proteins & nuclei segments 
 c) Resized Combined Lamin-B1 proteins & nuclei segments *) (*refer to Supplementary material & Appendix A*)

(* Input cell 54 *)
{allProteins[blueimg_, redimg_]:=allProteins[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], redimg}, {"MaskedImage", "Centroid"}];, , LB1prolist=allProteins[laminB1Blue, laminB1Red];, , LB1prolistV=Values[LB1prolist];, , LB1proPics=LB1prolistV[[All, 1]];, , LB1rProtein=resnetF/@LB1proPics;, , Dimensions[LB1rProtein], }

(* Input cell 55 *)
{{LB1myPlot1a, LB1clus1a}=ProteinClusterAndPlot[LB1rProtein, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {LB1myPlot1b, LB1clus1b}=ProteinClusterAndPlot[LB1rProtein, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];, }

(* Input cell 56 *)
Grid[{{Style["<Overlay with 2 Clusters of Lamin-B1 Proteins Only 
(KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Lamin-B1 Proteins Only 
(KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[laminB1Blue, laminB1Red, laminB1Grey, LB1clus1a], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[laminB1Blue, laminB1Red, laminB1Grey, LB1clus1b], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 57 *)
{allNuclei[blueimg_, redimg_]:=allNuclei[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], blueimg}, {"MaskedImage", "Centroid"}];, , LB1nucList=allNuclei[laminB1Blue, laminB1Red];, , LB1nuPics = Values[LB1nucList][[All, 1]];}

(* Input cell 58 *)
{LB1allNuPro=Transpose[{LB1nuPics, LB1proPics}];, , LB1rNuPro=resnetF/@LB1allNuPro;, , Dimensions[LB1rNuPro]}

(* Input cell 59 *)
{LB1mergedFeatures=Map[Join[# [[1]], # [[2]]]&, (*take nucleus (first) then protein (second)*) 
 LB1rNuPro];, , Dimensions[LB1mergedFeatures]}

(* Input cell 60 *)
{{LB1myPlot2a, LB1clus2a}=NucleusProteinClusterAndPlot[LB1mergedFeatures, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {LB1myPlot2b, LB1clus2b}=NucleusProteinClusterAndPlot[LB1mergedFeatures, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];}

(* Input cell 61 *)
Grid[{{Style["<Overlay with 2 Clusters of Lamin-B1 Nuclei & Proteins 
(KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Lamin-B1 Nuclei & Proteins 
(KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[laminB1Blue, laminB1Red, laminB1Grey, LB1clus2a], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[laminB1Blue, laminB1Red, laminB1Grey, LB1clus2b], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 62 *)
{ResizeImageForResNet[img_]:=ColorConvert[ImageResize[RemoveAlphaChannel[img], {224, 224}], "RGB"];, , LB1nucleusResized=ResizeImageForResNet/@LB1nuPics;, , LB1proteinResized=ResizeImageForResNet/@LB1proPics;, }

(* Input cell 63 *)
{ExtractFeatures[img_]:=resnetF[img];, , LB1nucleusFeatures=ExtractFeatures/@LB1nucleusResized;, , LB1proteinFeatures=ExtractFeatures/@LB1proteinResized;, , LB1featureVectors=MapThread[Join, {LB1nucleusFeatures, LB1proteinFeatures}];, , Dimensions[LB1featureVectors]}

(* Input cell 64 *)
{{LB1myNPlot, LB1NproteinClusters}=NucleusProteinClusterAndPlot[LB1featureVectors, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {LB1myNPlot2, LB1NproteinClusters2}=NucleusProteinClusterAndPlot[LB1featureVectors, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];}

(* Input cell 65 *)
Grid[{{Style["<Overlay with 2 Clusters of Resized Lamin-B1 Nuclei & 
Proteins (KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Resized Lamin-B1 Nuclei & 
Proteins (KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[laminB1Blue, laminB1Red, laminB1Grey, LB1NproteinClusters], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[laminB1Blue, laminB1Red, laminB1Grey, LB1NproteinClusters2], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 66 *)
(*Replicating steps for spatial distribution above on 
 a) Nesprin-1 proteins only 
 b) Combined Nesprin-1 proteins & nuclei segments 
 c) Resized Combined Nesprin-11 proteins & nuclei segments *) 
 (*refer to Supplementary material & Appendix B*)

(* Input cell 67 *)
{allProteins[blueimg_, redimg_]:=allProteins[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], redimg}, {"MaskedImage", "Centroid"}];, , N1prolist=allProteins[nesprin1Blue, nesprin1Red];, , N1prolistV=Values[N1prolist];, , N1proPics=N1prolistV[[All, 1]];, , N1rProtein=resnetF/@N1proPics;, , Dimensions[N1rProtein], }

(* Input cell 68 *)
{{N1myPlot1a, N1clus1a}=ProteinClusterAndPlot[N1rProtein, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {N1myPlot1b, N1clus1b}=ProteinClusterAndPlot[N1rProtein, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];, }

(* Input cell 69 *)
Grid[{{Style["<Overlay with 2 Clusters of Nesprin-1 Proteins Only 
(KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Nesprin-1 Proteins Only 
(KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[nesprin1Blue, nesprin1Red, nesprin1Grey, N1clus1a], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[nesprin1Blue, nesprin1Red, nesprin1Grey, N1clus1b], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 70 *)
{allNuclei[blueimg_, redimg_]:=allNuclei[blueimg, redimg]=ComponentMeasurements[{Erosion[Dilation[MorphologicalComponents[segmentation[blueimg]], 10], 5], blueimg}, {"MaskedImage", "Centroid"}];, , N1nucList=allNuclei[nesprin1Blue, nesprin1Red];, , N1nuPics = Values[N1nucList][[All, 1]];}

(* Input cell 71 *)
{N1allNuPro=Transpose[{N1nuPics, N1proPics}];, , N1rNuPro=resnetF/@N1allNuPro;, , Dimensions[N1rNuPro]}

(* Input cell 72 *)
{N1mergedFeatures=Map[Join[# [[1]], # [[2]]]&, (*take nucleus (first) then protein (second)*) 
 N1rNuPro];, , Dimensions[N1mergedFeatures]}

(* Input cell 73 *)
{{N1myPlot2a, N1clus2a}=NucleusProteinClusterAndPlot[N1mergedFeatures, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {N1myPlot2b, N1clus2b}=NucleusProteinClusterAndPlot[N1mergedFeatures, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];}

(* Input cell 74 *)
Grid[{{Style["<Overlay with 2 Clusters of Nesprin-1 Nuclei & Proteins 
(KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Nesprin-1 Nuclei & Proteins 
(KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[nesprin1Blue, nesprin1Red, nesprin1Grey, N1clus2a], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[nesprin1Blue, nesprin1Red, nesprin1Grey, N1clus2b], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]

(* Input cell 75 *)
{ResizeImageForResNet[img_]:=ColorConvert[ImageResize[RemoveAlphaChannel[img], {224, 224}], "RGB"];, , N1nucleusResized=ResizeImageForResNet/@N1nuPics;, , N1proteinResized=ResizeImageForResNet/@N1proPics;, }

(* Input cell 76 *)
{ExtractFeatures[img_]:=resnetF[img];, , N1nucleusFeatures=ExtractFeatures/@N1nucleusResized;, , N1proteinFeatures=ExtractFeatures/@N1proteinResized;, , N1featureVectors=MapThread[Join, {N1nucleusFeatures, N1proteinFeatures}];, , Dimensions[N1featureVectors]}

(* Input cell 77 *)
{{N1myNPlot, N1NproteinClusters}=NucleusProteinClusterAndPlot[N1featureVectors, 2, Method->"KMeans", DistanceFunction->EuclideanDistance];, , {N1myNPlot2, N1NproteinClusters2}=NucleusProteinClusterAndPlot[N1featureVectors, 3, Method->"KMeans", DistanceFunction->EuclideanDistance];}

(* Input cell 78 *)
Grid[{{Style["<Overlay with 2 Clusters of Resized Nesprin-1 Nuclei & 
Proteins (KMeans, Level=1, Euclidean Distance)", Bold, 14], Style["<Overlay with 3 Clusters of Resized Nesprin-1 Nuclei & 
Proteins (KMeans, Level=1, Euclidean Distance)", Bold, 14]}, {Show[segmentationOverlayWithClustersUniversal[nesprin1Blue, nesprin1Red, nesprin1Grey, N1NproteinClusters], ImageSize->550], Show[segmentationOverlayWithClustersUniversal[nesprin1Blue, nesprin1Red, nesprin1Grey, N1NproteinClusters2], ImageSize->550]}}, Alignment->Center, Spacings->{2, 1}, Frame->All, FrameStyle->LightGray, ItemStyle->Directive[FontFamily->"Arial", FontSize->12]]