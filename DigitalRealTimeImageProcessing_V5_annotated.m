clear;
close all;

%Automated analyzer for real-time isothermal amplification
%by Erik Jue, 2018

%This MATLAB script performs the following functions
%Loads a 12-bit .tif image sequence for digital LAMP experiment
%Uses first image frame (low temp) to detect total # wells
%Uses 2nd to last image frame (end of experiment) to detect positive wells
%Track the intensity of the positive wells over time
%Apply gaussian smoothing and baseline subtraction
%Saves the data
%Repeats for each tiff stack in the folder 

%Requires Bio-Formats 5.8.1 by OME
%https://www.openmicroscopy.org/bio-formats/
%Requires Control ZEN Blue and the microscope from MATLAB (ReadImage6D.m) by Sebastian Rhode
%https://www.mathworks.com/matlabcentral/fileexchange/50079-control-zen-blue-and-the-microscope-from-matlab?focused=6875595&tab=function

%Notes specific to this data set
%Used a connectivity of 4 (cardinal directions only) for well detection
%Ignore last frame since it does not capture the full exposure
%Sample fluorescence is high at low temp, use this to detect all wells
%Sample fluorescence decreases with high temp, ignore all heating frames 

%% VARIABLES FOR EDITING

%FOLDER TO ANALYZE, WILL ANALYZE ALL .TIF FILES IT FINDS
folder = 'C:\Users\Lab User\Desktop\test'; %'C:\user\data';

%SAVE LOCATION
saveHere = 'C:\Users\Lab User\Desktop\Results'; %'C:\user\results\';

%Adds the trailing \ if applicable
if saveHere(end)~='\'
    saveHere = strcat(saveHere, '\');
end

%SETTINGS, 1 means on, 0 means off
%Set to 1 to save images
imageSave = 1;
%Set to 1 to save summary info to excel file
excelSaveSummary = 1;
%Set to 1 to save intensity curves to excel file
excelSaveIntensity = 1;

%baselineAdjust averages the intensities between all frames in between 
%baselineStart and baselineEnd and subtracts this value for all frames
%Set to 1 to enable background subtraction
baselineAdjust = 1;

%This is the frame number when the sample reaches the correct temp
frameAtTemp = 6;
baselineStart = frameAtTemp;
baselineEnd = frameAtTemp+5;

%slopeAdjust finds the slope between the baselineEnd frame and
%baselineStart frame and propagates the slope correction for all frames
%This can be used to account for slow drift of negatives wells
%Note: do not use if raw data is very noisy
%Set to 1 to enable background slope correction
slopeAdjust = 0;

%Well Detection Settings
%Multi-level thresholding
%This techniques applies a threshold to detect a set of wells and selects
%selects wells that fall in the areaBound and majorAxisBound. Multiple
%rounds of different thresholds are applied to detect a greater number of
%wells

%Manually determined threshold levels for first image
im1_thresh = [.07 .08 .09 .1 .11 .12 .13 .14 .15 .175 .2 0.3 0.35 0.4 0.5];

%Manually determine threshold levels for the masking image
mask_thresh = [.05 .075 .1 .125 .15 .175 .2 .25 .3 .4];

%lower and upper cutoff for the area of each well
%Units in # of pixels
areaBound = [4 12];

%lower and upper cutoff for the major axis of each well (ensures circularity)
majorAxisBound = [2 5];

%Thresholding is used for automated time-to-positive determination
%Will calculate ttp after baselineAdjust and slopeAdjust if applicable
threshold = 100;

%A gaussian filter is used to smooth the intensity curves
%Set the gaussian window smoothing size
gaussWinSize = 10;

%Maximum allowable slope to detect
maxSlope = 200;

%Threshold for the max slope required to call a well positive
maxSlopeThreshold = 30;

%% PROGRAM START HERE

%Identify all the .tif files
myFiles = dir(fullfile(folder, '*.tif'));

%For each .tif, run image analysis
for a=1:size(myFiles,1)
    %Generate the specific file name
    baseFileName = myFiles(a).name;
    filename = fullfile(folder, baseFileName);
    
    %Load the image into memory
    tifStack=ReadImage6D(filename);
    disp(['Loaded ', filename])
    
    %The data file tifStack stores image info in {1} and metadata in {2}
    %Grab the metadata
    numImages = tifStack{2}.SizeZ;
    %For our images, %Y pixels are actually X pixels
    xpixels = tifStack{2}.SizeY; 
    ypixels = tifStack{2}.SizeX;
    
    %tifStack stores image info {1} in the format below
      %1 = Series
      %2 = SizeC
      %3 = SizeZ
      %4 = SizeT
      %5 = SizeX
      %6 = SizeY
    %Depending on how the .tif was generated, image # may be in slot 3 or 4

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %COUNT TOTAL WELLS USING 1ST IMAGE
    %Matlab does not natively support 12-bit, convert to 16-bit
    im1 = cast(16*squeeze(tifStack{1}(1,1,1,1,:,:)), 'uint16');
    
    %Detect all wells for the first image
    im1_bestMask = zeros(size(im1));
    for k = 1:size(im1_thresh,2)
        bwthresh = imbinarize(im1, im1_thresh(k));

        %Apply area filter
        areaFilter = bwpropfilt(bwthresh, 'Area', areaBound, 4);
        %Apply major axis filter
        majorAxisFilter = bwpropfilt(areaFilter, 'MajorAxisLength', majorAxisBound, 4);
        
        im1_bestMask = im1_bestMask + majorAxisFilter;
    end
    im1_bestMask = im1_bestMask~=0;
    
    %Calculate the total number of wells using the mask
    im1_result = bwconncomp(im1_bestMask, 4);
    totalNumWells = im1_result.NumObjects;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %COUNT POSITIVE WELLS USING 2ND TO LAST IMAGE
    maskIm = cast(16*squeeze(tifStack{1}(1,1,numImages-1,1,:,:)), 'uint16');
    
    bestMask = zeros(size(im1));
    for k = 1:size(mask_thresh, 2)
        bwthresh = imbinarize(maskIm, mask_thresh(k));

        %Apply area filter
        areaFilter = bwpropfilt(bwthresh, 'Area', areaBound, 4);
        %Apply major axis filter
        majorAxisFilter = bwpropfilt(areaFilter, 'MajorAxisLength', majorAxisBound, 4);
        
        bestMask = bestMask + majorAxisFilter;
    end
    bestMask = bestMask~=0;

    %Generate the mask to use for all images
    labeledmask= bwconncomp(bestMask, 4);
    disp('Mask found, begin avg intensity calculations')

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FIND INTENSITY OF ALL POSITIVE WELLS OVER TIME
    
    %Call regionprops to count the number of positive wells
    resultTemp = regionprops(labeledmask, bestMask, 'MeanIntensity');
    numWells = size(resultTemp, 1);
    
    %Initialize matrix to track all intensity curves
    intensity = zeros(numImages-1, numWells);
    
    %Loop through all images and find the avg intensity of each well
    for i = 1:numImages-1
        if mod(i,10)==0
            disp(strcat('Calculating image #', num2str(i)))
        end

        current_im_u16 = cast(squeeze(tifStack{1}(1,1,i,1,:,:)), 'uint16');

        result = regionprops(labeledmask, current_im_u16, 'MeanIntensity', 'Area');
        intensity(i, :) = cat(1, result.MeanIntensity);
    end
    
    %Pull out the area for each well
    wellArea = zeros(numWells, 1);
    wellArea(:) = cat(1, result.Area);
    disp('All well intensities found')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DATA SMOOTHING, BKGND CORRECTION, AND SLOPE CORRECTION
    %Apply heuristically determined moving average window filter
    smooth = smoothdata(intensity, 'gaussian', gaussWinSize);

    %Background average subtraction
    baseline = zeros(1, numWells);
    baselineData = intensity;

    if baselineAdjust==1
        %Average the intensities between baselineStart and baselineEnd for
        %each well
        for i = 1:numWells
            baseline(i) = mean(smooth(baselineStart:baselineEnd, i));
        end
        
        %Subtract the baseline for each well
        for i = 1:numWells
            for j = 1:numImages-1
               baselineData(j,i) = smooth(j,i) - baseline(i); 
            end
        end
    end
    
    %Adjust the slope from baseline
    if slopeAdjust==1
        for i = 1:size(baselineData,1)
            slope = (baselineData(baselineEnd, i) - baselineData(baselineStart, i)) / (baselineEnd - baselineStart + 1);
            baselineSlope = linspace(0, (numImages-1)*slope, numImages);
            baselineSlope = baselineSlope - (((baselineEnd+baselineStart)/2)-1)*slope;

            baselineData(i, :) = baselineData(i, :) - baselineSlope(i);
        end
    end
    disp('Finished data smoothing')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DATA PROCESSING
    
    %Determine which wells cross the intensity threshold
    %Initialize an array to hold time-to-positive values for each well
    ttp_thresh = zeros(numWells,1 );

    %Do not count any time-to-positives before sample has reached temp
    for i=frameAtTemp:numImages-1
        for j=1:numWells
            %Record the first frame that a well crosses the threshold
            if (baselineData(i,j)>threshold && ttp_thresh(j)==0)
                ttp_thresh(j) = i;
            end
        end
    end
    
    %Generate new array that only contains wells that cross the threshold
    ttp_thresh_clean = ttp_thresh(ttp_thresh>0);
    
    %Generate array of number of wells on per frame
    Intensity_Histogram = zeros(numImages,1);
    for i=1:numImages
        Intensity_Histogram(i) = sum(ttp_thresh_clean(:)==i);
    end
    
    %Find the max slope for each well
    fast_Slope = zeros(numWells, 1);
    for p = 1:numWells
        fastestSlope = max(diff(baselineData(:,p)));
        %Ignores extremely high slopes that are imaging artifacts
        if fastestSlope < maxSlope
            fast_Slope(p) = fastestSlope;
        end
    end
    
    %Determine which wells cross the intensity AND slope thresholds
    ttp_thresh_slope = ttp_thresh(fast_Slope>maxSlopeThreshold);
    ttp_thresh_slope_clean = ttp_thresh_slope(ttp_thresh_slope>0);
    
    Intensity_Slope_Histogram = zeros(numImages, 1);
    for i=1:numImages
        Intensity_Slope_Histogram(i) = sum(ttp_thresh_slope_clean(:)==i);
    end
    
    %Determine which wells cross slope threshold
    baselineDiffData = diff(baselineData,1,1);
        
    ttp_slopeOnly = zeros(numWells,1);
    for q=frameAtTemp:numImages-2
        for r=1:numWells
            %Record the first frame that a well crosses the threshold
            if (baselineDiffData(q,r)>maxSlopeThreshold && ttp_slopeOnly(r)==0)
                ttp_slopeOnly(r) = q;
            end
        end
    end
    
    %Generate new array that only contains wells that cross the threshold
    ttp_slopeOnly_clean = ttp_slopeOnly(ttp_slopeOnly>0);
    
    Slope_Histogram = zeros(numImages-1, 1);
    for i=1:numImages-1
        Slope_Histogram(i) = sum(ttp_slopeOnly_clean(:)==i);
    end
    
    %Find the max intensity of each positive well
    endIntensity = zeros(numWells, 1);
    for m = 1:numWells
        endIntensity(m) = (baselineData(numImages-1,m));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SAVE DATA
    %Break the filename into parts
    [filepath, name, ext] = fileparts(filename);
    
    %Create a folder to save the results
    if imageSave == 1
        disp('Saving Images')
        mkdir(strcat(saveHere, 'Results_', name));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot the generation of the masks
    figure(1)
    subplot(2,2,1);
    imshow(im1);
    title('First Image')
    
    subplot(2,2,2);
    imshow(im1_bestMask);
    title('Well Counting')
    
    subplot(2,2,3);
    imshow(maskIm);
    title('Masking Image')
    
    subplot(2,2,4);
    imshow(bestMask);
    title('Final Mask')
    
    if imageSave == 1
        saveas(gcf,strcat(saveHere, 'Results_', name,'\Thresholding'),'png')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Plot baselined intensities over time
    figure(2)
    
    %Downsample here to shown fewer curves
    plot(baselineData(:, :)); 
    
    %Set the plot axes
    axis([0 numImages -100 1000]);
    title([name ' Corrected Well Intensity']);
    xlabel('Frame Number (2 per minute)');
    ylabel('RFU');
    
    if imageSave == 1
        saveas(gcf,strcat(saveHere, 'Results_', name,'\BaselineIntensity'),'png')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Plot derivative of baseline over time
    figure(3)
    
    %Downsample here to shown fewer curves
    plot(baselineDiffData(:, :)); 
    
    %Set the plot axes
    axis([0 numImages -30 100]);
    title([name ' Corrected Derivative']);
    xlabel('Frame Number (2 per minute)');
    ylabel('Delta RFU');
    
    if imageSave == 1
        saveas(gcf,strcat(saveHere, 'Results_', name,'\BaselineDeriv'),'png')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%PLOT HISTOGRAM AND CDF FOR WELLS THAT CROSS INTENSITY THRESH
    
    %Only plot if the are wells that turned positive. 
    if ttp_thresh_clean > 1    
        figure(4)
        h=histcounts(ttp_thresh_clean, max(max(ttp_thresh_clean)) - min(min(ttp_thresh_clean)) + 1);

        %LineSpec Settings: line, marker, and colour
        plot(min(min(ttp_thresh_clean)):max(max(ttp_thresh_clean)), h, 'b');
        axis([0 numImages -inf inf]);
        title([name ' Intensity Histogram'])
        ylabel('Wells Turning On')
        xlabel('Frame Number (2 per Minute)')
        if imageSave == 1
            saveas(gcf,strcat(saveHere, 'Results_', name,'\Intensity_Hist'),'png')
        end
        
        %%Plot ttp CDF
        figure(5)
        intensity_cdfhandle = cdfplot(ttp_thresh_clean);
        axis([0 numImages 0 1]);
        title([name ' Intensity CDF, ', num2str(size(ttp_thresh_clean,1)), ' Wells'])
        xlabel('Frame Number (2 per Minute)')
        set(intensity_cdfhandle, 'LineStyle', '-', 'Color', 'r');
        
        if imageSave == 1
            saveas(gcf,strcat(saveHere, 'Results_', name,'\Intensity_CDF'),'png')
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%PLOT HISTOGRAM AND CDF FOR WELLS THAT CROSS INTENSITY AND SLOPE THRESH
    if ttp_thresh_slope_clean > 1    
        figure(6)
        h=histcounts(ttp_thresh_slope_clean, max(max(ttp_thresh_slope_clean)) - min(min(ttp_thresh_slope_clean)) + 1);

        %LineSpec Settings: line, marker, and colour
        plot(min(min(ttp_thresh_slope_clean)):max(max(ttp_thresh_slope_clean)), h, 'b');
        axis([0 numImages -inf inf]);
        title([name 'Intensity Slope Histogram'])
        ylabel('Wells Turning On')
        xlabel('Frame Number (2 per Minute)')
        if imageSave == 1
            saveas(gcf,strcat(saveHere, 'Results_', name,'\Intensity_Slope_Hist'),'png')
        end
        
        figure(7)
        intensity_slope_cdfhandle = cdfplot(ttp_thresh_slope_clean);  
        title([name ' Intensity Slope CDF, ', num2str(size(ttp_thresh_slope_clean,1)), ' Wells'])
        xlabel('Frame Number (2 per Minute)')
        set(intensity_slope_cdfhandle, 'LineStyle', '-', 'Color', 'r');

        if imageSave == 1
            saveas(gcf,strcat(saveHere, 'Results_', name,'\Intensity_Slope_CDF'),'png')
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%PLOT HISTOGRAM AND CDF FOR WELLS THAT SLOPE ONLY THRESH
    if ttp_slopeOnly_clean > 1    
        figure(8)
        h=histcounts(ttp_slopeOnly_clean, max(max(ttp_slopeOnly_clean)) - min(min(ttp_slopeOnly_clean)) + 1);

        %LineSpec Settings: line, marker, and colour
        plot(min(min(ttp_slopeOnly_clean)):max(max(ttp_slopeOnly_clean)), h, 'b');
        axis([0 numImages -inf inf]);
        title([name 'Intensity Slope Histogram'])
        ylabel('Wells Turning On')
        xlabel('Frame Number (2 per Minute)')
        if imageSave == 1
            saveas(gcf,strcat(saveHere, 'Results_', name,'\Slope_Hist'),'png')
        end
        
        figure(9)
        slope_cdfhandle = cdfplot(ttp_slopeOnly_clean);  
        title([name ' Intensity Slope CDF, ', num2str(size(ttp_slopeOnly_clean,1)), ' Wells'])
        xlabel('Frame Number (2 per Minute)')
        set(slope_cdfhandle, 'LineStyle', '-', 'Color', 'r');

        if imageSave == 1
            saveas(gcf,strcat(saveHere, 'Results_', name,'\Slope_CDF'),'png')
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SAVE TO EXCEL

    %SUMMARY FILE
    %TTP_hist: Contains the number of wells that turned on for each frame
    %fast_TTP_hist: Contains the number of wells that turned on after
    %   filtering by slope threshold
    %ttp: Lists time-to-positive for each tracked well that
    %   crosses the treshold
    %fast_ttp: Lists time-to-positive for each tracked well that
    %   crosses the treshold and slope threshold
    %maxSlope: Maximum slope/frame for each well
    %endIntensity: Final intensity (after baseline) for each well
    %wellArea: Area of each well

    %Total wells: Total wells detected on frame 1
    %Detected wells: number of wells detected on 2nd to last frame
    %Positive wells: number of detected wells that cross threshold
    %Fast slope pos wells: number of positive wells that cross slope threshold

    col_header={'Intensity_ttp','Intensity_hist', 'Intensity_slope_ttp', 'Intensity_slope_hist', 'slope_ttp', 'slope_hist','maxSlope', 'endIntensity', 'wellArea', '', ''};
    col1 = cell(numWells, 1); col1(1:size(ttp_thresh_clean,1),1) = num2cell(ttp_thresh_clean);
    col2 = cell(numWells, 1); col2(1:size(Intensity_Histogram,1),1) = num2cell(Intensity_Histogram);
    col3 = cell(numWells, 1); col3(1:size(ttp_thresh_slope_clean,1),1) = num2cell(ttp_thresh_slope_clean);
    col4 = cell(numWells, 1); col4(1:size(Intensity_Slope_Histogram,1),1) = num2cell(Intensity_Slope_Histogram);
    col5 = cell(numWells, 1); col5(1:size(ttp_slopeOnly_clean,1),1) = num2cell(ttp_slopeOnly_clean);
    col6 = cell(numWells, 1); col6(1:size(Slope_Histogram,1),1) = num2cell(Slope_Histogram);
    
    col7 = cell(numWells, 1); col7(1:size(fast_Slope,1),1) = num2cell(fast_Slope);
    col8 = cell(numWells, 1); col8(1:size(endIntensity,1),1) = num2cell(endIntensity);
    col9 = cell(numWells, 1); col9(1:size(wellArea,1),1) = num2cell(wellArea);
    col10= cell(numWells, 1); col10(1:18,1) = {'Total wells'; 'Detected wells'; 'Intensity wells'; 'Intensity slope wells'; 'slope wells';...
        ''; 'Settings'; 'frameAtTemp'; 'baselineStart'; 'baselineEnd'; 'areaLowBound'; 'areaHighBound';...
        'majorAxisLowBound'; 'majorAxisHighBound'; 'threshold'; 'GaussWindow'; 'maxSlope'; 'maxSlopeThreshold'};
    col11 = cell(numWells, 1); col11(1:18,1) = num2cell([totalNumWells; numWells; size(ttp_thresh_clean,1); size(ttp_thresh_slope_clean,1); size(ttp_slopeOnly_clean,1);...
        0; 0;          frameAtTemp; baselineStart; baselineEnd; areaBound(1); areaBound(2);
        majorAxisBound(1); majorAxisBound(2); threshold; gaussWinSize; maxSlope; maxSlopeThreshold]);
    col11(6:7,1)={''; ''};

    saveSummary = [col_header; col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11];

    if excelSaveSummary == 1
        disp('Saving Summary Excel File')
        xlswrite(strcat(saveHere, 'Results_', name,'\',name,' Summary.xlsx'), saveSummary)
    end 

    %INTENSITY FILE
    %Total number of wells stored as the name of the sheet in Intensity
    %Intensity: raw traces for each tracked well
    %baselineData: traces after gaussian smoothing, baseline avg
    %subtraction and slope correction (if applicable).
    if excelSaveIntensity == 1
        disp('Saving Intensity Excel File')
        xlswrite(strcat(saveHere, 'Results_', name,'\',name,' Intensity.xlsx'), transpose(intensity), strcat('Intensity_', num2str(totalNumWells)))
        xlswrite(strcat(saveHere, 'Results_', name,'\',name,' Intensity.xlsx'), transpose(baselineData), 'BaselineCorrected')
        xlswrite(strcat(saveHere, 'Results_', name,'\',name,' Intensity.xlsx'), transpose(baselineDiffData), 'BaselineDerivative')
    end 

    disp('Done Processing')
end