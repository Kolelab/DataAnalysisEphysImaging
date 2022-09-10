%% DaMatlab.m
% This script takes electrophysiological and optically recorded data (using 
% the neuroplex software) and analyses it together. You can select a ROI 
% based on a .det file (selected in Neuroplex). The two dataset are
% temporally aligned (based on the threshold in the ephys data). The the
% light periods are selected, the trials are averaged and a corresponding
% time vector is also created (with time of threshold t = 0 ms). Then the
% optical data is collapsed over x and corrected for bleach. Data is
% exported and further analyzed in ExtractParams.m.

% Organize data in a folder named date_cell_condition, with condition
% described in 4 characters. Example: 220212_c1_ctrl
% In that folder have the .da files (optical data with Neuroplex) for all
% trials seperate. The ephys trials are combined in Axograph and then
% exported as a matlab file, also store that in this folder. Then in the
% folder one level up have a folder called DetFiles where you store the
% .det files (containing the ROI pixels selected in Neuroplex) with the
% name: date_cell_roi.det

% At the end there are also two ways to plot the data. 


%% Select .da file

project = input('What is the project? 1 = AIS, 2 = Epi vs HG, 3 = Nod 4 = Calibration');
if project == 1
    project = 'AIS';
elseif project == 2
    project = 'HLG';
elseif project == 3
    project = 'Nod';
elseif project == 4
    project = 'Cal';
end


[file, pathfolder] = uigetfile(strcat('\\vs01\AXS_DATA\APshape\',project,'\tfa\*.da'),'Select a .da file','MultiSelect','on');
cellname = pathfolder(33:41);
condition = pathfolder(43:end-1);
mkdir(strcat('\\vs01\AXS_DATA\APshape\',project,'\BCData\',cellname, '_', condition))
outputpath = strcat('\\vs01\AXS_DATA\APshape\',project,'\BCData\',cellname, '_', condition,'\');

%% Select Axg file
% first combine episodes in axograph and then export these files to matlab

[Vm_RMP, dVm, Vm_APPeak, InjCur, Vm_th, Vm_ThresholdInd, VoltageSoma, Delay, InjToTHAxg, OffsetStimDelay, Sag] = analysisEphys(pathfolder, cellname, condition);

if length(file) == length(Delay)
    disp('Yay! Number of episodes in NP and Ephys match')
else
    disp('Oops! Mismatch in number of episodes between NP and Ephys')
end

disp('Ephys imported')

%% HArr         Read header-array
% The .da file is one long array that consists of several parts.
% The first part (the header-array or HArr) contains information about the recording, such as the  dimensions.
% HArr has the fixed size of (1,2560).
% We read that part of the file first to determine these parameters. 

if iscell(file)>0.5
    filename = sprintf(file{1,1});
else
    filename = sprintf(file);
end

fileID = fopen(fullfile(pathfolder,filename),'r');
HArr = fread(fileID, 2560, 'integer*2');
fclose(fileID);

clear fileID
clear ans
clear filename

% After HArr comes FuAr (functional array), whose size depends on the
% amount of data. So we look at HArr to find the dimensions of the data:

x = HArr(386);
y = HArr(385);
z = HArr(5);
FrameRate = 1/(HArr(389)*HArr(391)/1000000);      % HArr(389)*HArr(391) = duration of a frame in microseconds
EndFr = 2560+(x*y*z);                   
NrTrials = length(file);

disp('HARR imported')

%% Load ROI and BGR pixels
% ROI = all pixels that I selected in NP (will be analysed)
% BGR = background will be subtracted
% pixels are selected in Neuroplex and saved as a .det file. Load them in here.
% BGR is only necessary if you want dFF, I dont use it.

[ROIDetFile, ROIpath] = uigetfile(strcat('\\vs01\AXS_DATA\APshape\',project,'\tfa\DetFiles\', cellname,'_roi.det'),'Select a ROI .det file');
%[BGRDetFile, BGRpath] = uigetfile(strcat('\\vs01\AXS_DATA\APshape\',project,'\tfa\DetFiles\', cellname,'_bgr.det'),'Select a BGR .det file');

fileID = fopen(fullfile(ROIpath,ROIDetFile),'r');
ROIPixels = fscanf(fileID, '%f');
fclose(fileID);

ROIPixelNrs = reshape(1:x*y, [26 4]);
ROIPixelNrs(ROIPixels) = 0;
ROIPixerNrsT  = transpose(ROIPixelNrs);
[ROIPixelsX, ROIPixelsY] = find(ROIPixerNrsT == 0);    

% BGRPixelNrs = reshape(1:x*y, [26 4]);
% BGRPixelNrs(BGRPixels) = 0;
% BGRPixerNrsT  = transpose(BGRPixelNrs);
% [BGRPixelsX, BGRPixelsY] = find(BGRPixerNrsT == 0);  

clear ans fileID ROIpath  
clear ROIPixels ROIPixelNrs ROIPixerNrsT ROIPixerNrsTF ROIDetFile
% clear BGRPixels BGRPixelNrs BGRPixerNrsT BGRPixerNrsTF BGRDetFile
disp('ROI pixels imported')

%% FuAr     Reading Functional Array from multiple files 
% Now we can read the 'actual' data (which is in the second part; FuAr) and
% third part (BNC data, contains Ephys traces), all BNC's are stored 
% sequentially, HArr(5) contains number of frames

%BNC1 = FuAr3(1:HArr(5));                  
%BNC2 = FuAr3(HArr(5)+1:2*HArr(5));
%BNC3 = Ephys_temp(2*HArr(5)+1:3*HArr(5));
%BNC4 = FuAr3(3*HArr(5)+1:4*HArr(5));
%BNC5 = FuAr3(4*HArr(5)+1:5*HArr(5));
%BNC6 = FuAr3(5*HArr(5)+1:6*HArr(5));
%BNC7 = FuAr3(6*HArr(5)+1:7*HArr(5));
%BNC8 = FuAr3(7*HArr(5)+1:8*HArr(5));

FuAr = zeros(x,y,z,NrTrials);
NPEphys = zeros(z,NrTrials);

for file_i = 1:NrTrials
    filename = sprintf(file{1,file_i});
    fileID = fopen(fullfile(pathfolder,filename),'r');
    Ephys_temp = fread(fileID, Inf, 'integer*2');
    Ephys_temp(1:EndFr) = [];                                 % remove HArr and FuAr
    NPEphys(:,file_i) = Ephys_temp(HArr(5)+1:2*HArr(5));      % change this to correct BNC
    clear Ephys_temp;
   
    fileID = fopen(fullfile(pathfolder,filename),'r');
    FuArTemp = fread(fileID, EndFr, 'integer*2');
    FuArTemp(1:2560) = [];                                    % Remove HArr from FuAr
    FuArTemp = reshape(FuArTemp, z, y, x);                    % Data is stored sequentially per pixel
    FuArMatr = permute(FuArTemp, [3, 2, 1]);                  % Reshape and permute the data into a matrix, with time in the 3rd dimension
    FuAr(:,:,:,file_i) = FuArMatr;
    fclose(fileID);
    clear FuArTemp FuArMatr fileID;
end

clear file_i ans filename
disp('FuAR imported')

%% FuArAligned      Temporal alignment to Ephys traces
% we need to temporally align all trials: use delay from Ephys

DelayR = round(Delay/5)+1;
DelayEnd = max(DelayR) - (DelayR);

FuArAligned = nan(size(FuAr,1), size(FuAr,2),size(FuAr,3)-max(DelayR), size(FuAr,4));
NPEphys = NPEphys';
EphysAligned = nan(size(NPEphys,1), size(NPEphys,2)-max(DelayR));

for trial = 1:NrTrials
    EphysTemp = NPEphys(trial,:);
    EphysAligned(trial,:) = EphysTemp(DelayR(trial):end - DelayEnd(trial) - 1);
    FuArTemp = FuAr(:,:,:,trial);
    FuArAligned(:,:,:,trial) = FuArTemp(:,:,DelayR(trial):end - DelayEnd(trial) - 1);
end

NPEphysMean = mean(EphysAligned,1);

% Figure if you want to check if alignment is correct
% figure
% hold on
% plot(1:1:length(NPEphys),NPEphys,'k')
% plot(1:1:length(EphysAligned),EphysAligned,'b')
% plot(1:1:length(EphysAligned),NPEphysMean,'r')
% xlim([1900 2200])   

disp('Data is temporally aligned')

%% FuArLight        Finding frame on and off sets 
% Light on and off are found as max and min in differentiation, but trace
% needs one round of filtering to remove noise which is amplified in diff

% Binomial filter 
BinomialCoeff = conv([1/2 1/2],[1/2 1/2]); 
FuArFilt = filter(BinomialCoeff, 1, FuArAligned,[],3);      % filter with 1 pass binomial along 3rd dimension
FuArFilt(:,:,1,:) = [];                                     % filter gives a 'delay' correct for it by removing first datapoint
FuArFilt(:,:,1,:) = FuArAligned(:,:,1,:);                   % make first and last the same as unfiltered
FuArFilt(:,:,end+1,:) = FuArAligned(:,:,end,:);
diffFuArFilt = diff(FuArFilt,1,3);

% find the two/three highest and lowest peaks for start and end of light periods 
% On should be from periods with most delay and vice versa, so select
% shortest light period
[~, leastDelayedTrial] = min(DelayR);
temp = nan(9999,1);
temp = squeeze(mean(mean(diffFuArFilt(ROIPixelsX,ROIPixelsY,:,leastDelayedTrial),1,'omitnan'),2,'omitnan'));
[~,OnLocs,~,~] = findpeaks(temp,'SortStr','descend'); 
OnLocs(3:end) = [];
OnLocs = sort(OnLocs);

[~, mostDelayedTrial] = max(DelayR);
temp = squeeze(mean(mean(diffFuArFilt(ROIPixelsX,ROIPixelsY,:,mostDelayedTrial),1,'omitnan'),2,'omitnan'));
[~,~,~,PeakProm] = findpeaks(temp*-1); 
[~,OffLocs,~,~] = findpeaks(temp*-1,'SortStr','descend'); 
OffLocs(3:end) = [];
OffLocs = sort(OffLocs);

% figure
% subplot(3,1,1);
% hold on
% plot(squeeze(FuArAligned(ROIPixelsX(1),ROIPixelsY(1),:,:)))
% plot(OnLocs(1),0,'ro')
% plot(OffLocs(1),0,'bo')
% xlim([OnLocs(1)-100 OffLocs(1)+100])
% subplot(3,1,2);
% hold on
% plot(squeeze(FuArAligned(ROIPixelsX(1),ROIPixelsY(1),:,:)))
% plot(OnLocs(2),0,'ro')
% plot(OffLocs(2),0,'bo')
% xlim([OnLocs(2)-100 OffLocs(2)+100])
% subplot(3,1,3);
% plot(VoltageSoma)
% xlim([OnLocs(1)*5-100 OffLocs(1)*5+100])

% I do the AP in the first (or second) light period and the calibration in 
% second (or third) and take five frames for rise and fall of light
OnLoc = OnLocs(1)+5;        % change depending on AP in first or second light period
OffLoc =  OffLocs(1)-5;
CalOnLoc = OnLocs(2)+5;
CalOffLoc = OffLocs(2)-5;

FuArLight = FuArAligned(:,:,OnLoc:OffLoc,:);
FuArHyper = FuArAligned(:,:,CalOnLoc:CalOffLoc,:);     % update this based on ephys protocol
NrLightFrames = length(FuArLight);
NrLightFramesCal = size(FuArHyper,3);
NPEphysLight = NPEphysMean(OnLoc:OffLoc);

FrameNumber = 0:1:NrLightFrames-1;
FrameTime = FrameNumber / (FrameRate / 1000);
FrameNumberCal = 0:1:NrLightFramesCal-1;
FrameTimeCal = FrameNumberCal / (FrameRate / 1000);

disp('Light ON periods are selected')

%% FuArMean     Average over trials

FuArMean = mean(FuArLight, 4);
FuArCalMean = mean(FuArHyper, 4);

disp('Trials are averaged')

%% Define Time for Ephys and imaging data

% Interpolate the NPEphys, so the traces can be aligned to the peak better
NPEphysMean100Hz = spline(1:1:length(NPEphysMean),NPEphysMean,1:0.2:length(NPEphysMean));
[~, NPVm_APPeakInd] = max(NPEphysMean100Hz);
[~, Vm_APPeakInd] = max(VoltageSoma);

% Correct for small difference in time
Vm = VoltageSoma((OnLoc*5)-(NPVm_APPeakInd-Vm_APPeakInd):(OffLoc*5)-(NPVm_APPeakInd-Vm_APPeakInd));

% with these lines I checked that the alignment is perfect
% NPE = NPEphysMean100Hz(OnLoc*5:OffLoc*5); 
% a = [rescale(Vm')  rescale(NPE')]; 
% plot(a)

% Then find the time to correspond to the voltage data
VmTime = (0:1:length(Vm)-1)/100;
Vm_dVdt = diff(Vm) ./ diff(VmTime);

% Set Vm_Threshold time = 0
Vm_ThresholdInd = find(Vm_dVdt > 50, 1,'first');
VmTime = VmTime - (Vm_ThresholdInd/100)+0.01;
FTime = (0:1:length(FuArMean)-1)/20;
FTime = FTime - (Vm_ThresholdInd/100)-0.02;
FrameTimeCalUps = FrameTimeCal(1):0.01:FrameTimeCal(end);

% here you can check alignment:
% hold on
% plot(FTime, rescale(NPEphysLight))
% plot(VmTime,rescale(Vm))

% find onset stim. Vm_ThresholdInd is the first, 
% so this is the average onsetstim, will be off with temporal jitter
OnsetStimInd = Vm_ThresholdInd - InjToTHAxg;
OnsetStim = VmTime(Vm_ThresholdInd - InjToTHAxg);
FOnsetStim = find(FTime < OnsetStim,1, 'last');

% Assume the AP occurs max 1 ms before or 2 ms after somatic AP
APLimMin = Vm_ThresholdInd-100;
APLimMax = Vm_ThresholdInd+200;

clear FrameTime FuArCalFilt FuArCalMeanROI FuArCalMeany FuArCalMeanyPDI 
clear FuArFilt FuArMeanROI FuArMeany FuArFilt FuARMean FuArMeanyPDI FuArTemp FrameRate
clear OnsetStimAxg OffsetStimDelay

disp('VmTime and VFTime created')

%% Collapse over x
% Optional. I collapse over x because that works with the axon.

FuArMeanROI = nan(size(FuArMean));
FuArCalMeanROI = nan(size(FuArCalMean));

for i_pix = 1:length(ROIPixelsX)
    FuArMeanROI(ROIPixelsX(i_pix),ROIPixelsY(i_pix),:) = FuArMean(ROIPixelsX(i_pix),ROIPixelsY(i_pix),:);
    FuArCalMeanROI(ROIPixelsX(i_pix),ROIPixelsY(i_pix),:) = FuArCalMean(ROIPixelsX(i_pix),ROIPixelsY(i_pix),:);
end

FuArMeany = squeeze(mean(FuArMeanROI, 1,'omitnan'));
FuArCalMeany = squeeze(mean(FuArCalMeanROI,1,'omitnan'));

clear FuArFilt FuArAligned EphysLight FrameNumber FuArLight
clear OnLocs OffLocs diffFuArFilt leastDelayedTrial mostDelayedTrial temp
clear EphysAligned CalOffLoc CalOnLoc Delay DelayEnd EphysTemp FuAr 
clear FuArCalMean EndFr FuArHyper FuArLight HArr NPEphys NrLightFrames
clear NrLightFramesCal PeakProm project ROIPixelsX ROIPixelsY Vm_ThresholdInd trial

disp('Data is 2D')

%% Bleach Correction calibration 
% BleachFit contains the fit parameters 

% first fit calibration

FuArCalMeanMean = mean(FuArCalMeany,1,'omitnan')';
FitTo = FuArCalMeanMean;
%FitTo(1:30) = nan;
% FitTo(end-20:end) = nan;
MeanFitTo = mean(FitTo,'omitnan');
xfit = (1:1:length(FitTo))';
idxValid = ~isnan(FitTo);
ft = fittype('a * exp(-x / t) + c'); 
fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Start',[10, MeanFitTo, 1],...  
    'Lower',[10, mean(FitTo,'omitnan')-500, 0],...
    'Upper',[200, mean(FitTo,'omitnan')-10, 200]);

BleachFitCal = fit(xfit(idxValid), FitTo(idxValid), ft, fo);
BleachFittedCal = feval(BleachFitCal,xfit);

hold on
plot(FuArCalMeanMean)
plot(FitTo)
plot(BleachFittedCal)


BleachFittedCalNorm = BleachFittedCal ./ BleachFittedCal(1,:);
BleachFittedCalNormM = permute(BleachFittedCalNorm,[2 1]);
BleachFittedCalNormM = repmat(BleachFittedCalNormM,y,1);

FuArCalBC = FuArCalMeany ./ BleachFittedCalNormM;

disp('Calibration pulse is bleach corrected')

%% Fit bleach correction AP

BlFitTo = 1:FOnsetStim-0; % generally 1:onsetstim. With a lot of temporal jitter in APs, onsetstim is not correct. Important to correct here, since the same range is used for calibration.
FuArMeanMean = mean(FuArMeany(1:y,:),1,'omitnan');
FitTo = nan(size(FuArMeanMean'));
FitTo(BlFitTo) = FuArMeanMean(BlFitTo);           
idxValid = ~isnan(FitTo);
xfit = (1:1:length(FitTo))';

ftbl = fittype('a * exp(-x / t) + c');                           % a+cis where is starts, c is where fit ends, t is steepness with low number is steeper
fobl = fitoptions('Method', 'NonlinearLeastSquares',...        
    'Start',[BleachFitCal.a, BleachFitCal.c, BleachFitCal.t],...
    'Lower',[0,   mean(FitTo,'omitnan')-150, 0],...          
    'Upper',[100, mean(FitTo,'omitnan')-0,   300]);              % sometimes i tinker with the Upperbound for c, so force the fit below the trace. Note in notes (a few lines down) what I adjusted.

                        
BleachFit = fit(xfit(idxValid), FitTo(idxValid), ftbl, fobl);
BleachFitted = feval(BleachFit,xfit);

figure
subplot(2,1,1)
hold on
plot(FTime, FuArMeanMean)
plot(FTime, FitTo)
plot(FTime, BleachFitted)
subplot(2,1,2)
hold on
plot(VmTime, Vm)
plot(OnsetStim, Vm(OnsetStimInd), 'x')

%% Bleach correction AP  

BleachFittedNorm = BleachFitted ./ BleachFitted(1,:);
BleachFittedNormM = permute(BleachFittedNorm,[2 1]);
BleachFittedNormM = repmat(BleachFittedNormM,y,1);

FuArBC = FuArMeany ./ BleachFittedNormM;


for i = 13:16
    figure
    subplot(2,1,1);
    hold on
    plot(FTime,squeeze(FuArMeany(i,:)),'b')
    plot(FTime,squeeze(FuArBC(i,:)),'r')
    plot(FTime(end) + FrameTimeCal,squeeze(FuArCalMeany (i,:)),'b')
    plot(FTime(end) + FrameTimeCal,squeeze(FuArCalBC(i,:)),'r')
    xlim([FTime(1)-1 (FTime(end) + FrameTimeCal(end) +1)])
    subplot(2,1,2);
    plot(VmTime, Vm)
    xlim([VmTime(1)-1 (VmTime(end) + FrameTimeCalUps(end) +1)])
    ylim([Vm_RMP-10 Vm_APPeak+10])
    pause
    close
end

Notes = cell(3,1);
Notes{1,1} = fobl;
Notes{2,1} = BlFitTo;
Notes{3,1} = input('What did you adjust in Bleach Fit?');

disp('Data is bleach corrected')



%% Export results

% Vm = Vm';
% VmTime = VmTime';
% FTime = FTime';
% FuArBC =FuArBC';
% CalF = CalF';


CalF = mean(FuArCalBC(:,1:end),2,'omitnan'); 
disp(strcat('You just saved:', ' ', cellname, '_', condition, '_all.mat'));
disp(' ');
save(strcat(outputpath,cellname, '_', condition, '_all.mat'))


clear




%% plot all pixel traces for 1 cell


colors = crameri('vik');
I = 1:256;
imagesc(I)
%colormap(colors)

figure
set(gcf, 'Position',  [100, 100, 450, 450])
for i = 1:26
    hold on
    plot(Time,VF_abs(i,:),'Color',colors(i*9,:))
end

plot([1 2], [0 0],'k')
plot([2 2], [0 50],'k')

plot(Time,Vm, 'k')
xlim([Time(APLimMin) Time(APLimMax)+3]) 

%saveas(gcf,strcat(outputpath,cellname, '_', condition,'_all example traces'),'epsc')

%% Space Time plot


%colors = cbrewer('div', 'RdYlBu', 256);

colors = crameri('vik');
%colors = flipud(colors);

ToPlot = [Vm; VF_abs]';
MaxP = max(max(ToPlot));
MinP = min(min(ToPlot));
ToPlotNorm = rescale(ToPlot);

hold on
subplot(1,2,1)
imagesc(ToPlotNorm(APLimMin:APLimMax,1),[0 1])
colormap(colors)
colorbar
yticklabels(Time(APLimMin):Time(APLimMax))
ylabel('time (ms)')
xlabel('soma')

daspect([1 10 1])


subplot(1,2,2)
imagesc(ToPlotNorm(APLimMin:APLimMax,2:end),[0 1])
colormap(colors)
caxis([-0.2 1])
xticklabels(17:17:85)
xlabel('distance from soma (\mum)')
yticklabels(Time(APLimMin:5:APLimMax))
ylabel('time (ms)')
daspect([1 10 1])
colorbar
%saveas(gcf,strcat(outputpath,cellname, '_', condition,'_spacetime'),'epsc')

% for i = 1:26
%     plot(FrameTime,dFF1Dnorm(:,i))
%     ylim([-1 1.5])
%     title(strcat('SNR = ',string(SNR1D(i))))
%     pause(0.5)
% end

end 







