%% Extract parameters from combined cells
% This script is used to extract parameters from patch-clamp and bleach 
% corrected imaging data. Bleach correct the data using DaMatlab.m
% This script loads the needed variables from the output of DaMatlab.m
% Requires Signal Processing Toolbox and Curve Fitting Toolbox

% general: all data corresponding to patch-clamp recordings is preceded by
% VM and all fluorescence data is preceded by F.
% This script was used to combine somatic patch-clamp recording and axonal
% fluorescence recordings. Regularly, the patch-clamp recording is referred
% to as 'soma' and the fluorescence recording is refered to as 'axon'.


%% Load Data
% This script required the output data from DaMatlab.m
% If your data is separated into projects, you can save them in separate
% folders. You need the following data structure:
% \\DataFolder\Project\BCdata
%                   ..\tfa

clear

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

% Mac
[BCDataPath]  = uigetfile_n_dir(strcat('/Volumes/AXS_DATA/APshape/',project,'/BCData/'),'select folders');

% VMware
%[BCDataPath]  = uigetfile_n_dir(strcat('\\VS01\AXS_DATA\APshape\',project,'\BCData\'),'select folders');

outputpath = '/Users/naomihanemaaijer/Dropbox/BK/';

NrExp = numel(BCDataPath);

cellId        = cell(NrExp,1);
condition     = cell(NrExp,1);
Vm            = cell(NrExp,1);
VmTime        = cell(NrExp,1);
FTime         = cell(NrExp,1);
FBCData       = cell(NrExp,1);
FCalInt       = cell(NrExp,1);
dVm           = cell(NrExp,1);
InjCur        = cell(NrExp,1);
NrTrials      = cell(NrExp,1);
PreStimFrames = cell(NrExp,1);
VmRMP         = cell(NrExp,1);
VmStartStimMS = cell(NrExp,1);

for iExp = 1:NrExp
    if project == 'Nod'
        cellId{iExp,1}        = BCDataPath{1,iExp}(end-15:end-7);
        condition{iExp,1}     = BCDataPath{1,iExp}(end-5:end);
    else
        cellId{iExp,1}        = BCDataPath{1,iExp}(end-13:end-5);
        condition{iExp,1}     = BCDataPath{1,iExp}(end-3:end);
    end
    filename              = strcat(BCDataPath{iExp}, '/', cellId{iExp}, '_', condition{iExp},'_all.mat');
    Vm{iExp,1}            = cell2mat(struct2cell( load(filename, 'Vm')));
    VmTime{iExp,1}        = cell2mat(struct2cell( load(filename, 'VmTime')));
    FTime{iExp,1}         = cell2mat(struct2cell( load(filename, 'FTime')));
    FBCData{iExp,1}       = cell2mat(struct2cell( load(filename, 'FuArBC')));
    FCalInt{iExp,1}       = cell2mat(struct2cell( load(filename, 'CalF')));
    dVm{iExp,1}           = cell2mat(struct2cell( load(filename, 'dVm')));
    InjCur{iExp,1}        = cell2mat(struct2cell( load(filename, 'InjCur')));
    NrTrials{iExp,1}      = cell2mat(struct2cell( load(filename, 'NrTrials')));
    PreStimFrames{iExp,1} = cell2mat(struct2cell( load(filename, 'BlFitTo')));
    VmRMP{iExp,1}         = cell2mat(struct2cell( load(filename, 'Vm_RMP')));
    VmStartStimMS{iExp,1} = cell2mat(struct2cell( load(filename, 'OnsetStim')));
end

y = 26;
SomaDist = (1.714:3.428:(y*3.428)-1.714)';

%% Optional: add averages of pixels at the end
% pixels to average, 4 pixels = 13.6 um, 3 pixels = 10.28, 5 = 17.14
% I tried different sizes of the spatial avering, since spatial averaging
% lowers the noise so increases the accuracy. However: you also increase
% variability because the AP is not constant over space. I compared 3, 4
% and 5 pixels. The SD of the FWHM in the axon was not very different 
% between these cases in dist AIS and internode , but in the prox AIS 
% smallest in the case of 4 pixels, so I chose that as the optimum between 
% spatial averaging and spatial accuracy. 

xas = ([1.714:(3.428):(y*(3.428))-1.714]);       % size of pixel = 3.428   
AvePixels = {3:6, 10:13, 19:22};               % average 4 pixels from 3 locations - proximal AIS, distal AIS and internode 
                                               % 1st = 6.856 - 20.568 um
                                               % 2nd = 30.852 - 44.56
                                               % 3rd = 61.704 - 75.416

                  
AveragedD = cell(NrExp,length(AvePixels));
AveragedC = cell(NrExp,length(AvePixels));
SASomaDist = nan(length(AvePixels), 1);

for iExp = 1:NrExp
    tempD = FBCData{iExp};
    tempC = FCalInt{iExp};
    MeanedD = nan(length(AvePixels), length(tempD));
    MeanedC = nan(length(AvePixels), size(tempC, 2));
    for iA = 1:length(AvePixels)
        MeanedD(iA,:)    = mean(tempD(AvePixels{iA},:),1,'omitnan');
        MeanedC(iA,:)    = mean(tempC(AvePixels{iA},:),1,'omitnan');
        SASomaDist(iA,:) = mean(SomaDist(AvePixels{iA},:),1,'omitnan');
    end
    FBCData{iExp} = [FBCData{iExp}; MeanedD];
    FCalInt{iExp} = [FCalInt{iExp}; MeanedC];
    clear MeanedD MeanedC tempD tempC
end

SASomaDist = [SomaDist; SASomaDist];
y = y + length(AvePixels);
xas  = [xas mean(xas(1,AvePixels{1},:),2,'omitnan') mean(xas(1,AvePixels{2},:),2,'omitnan') mean(xas(1,AvePixels{3},:),2,'omitnan')]; 

disp('Spatial averages added')

%% Calibrate AP signal 
% Calibrate F to Vm for every pixel
% FVm     = (Fbl ./ FdCal) .*dVm .* Decay; 
% FVm     = calibrated fluorescence difference from Vm_VmRMP as baseline
% FPreStim= mean F intensity before the stimulus
% Fbl     = F intensity with PreStim period as baseline = 0 
% FdCal   = difference in fluorescence in baseline and calibration pulse
% FCalInt = F intensity during cal pulse, loaded from DaMatlab.m
% dVm     = difference in mV from Vm_VmRMP to calibration pulse, loaded from DaMatlab.m
% Decay   = exponential decay of calibration pulse with distance
% Ouput is FVm (calibrated with RMP = 0) and VFmAbs (with RMP from soma)

if project == 'Nod'
    SomaDist = (1.714:3.428:(y*3.428)-1.714)' + 83;   %manually adjust the distance of the FOV in the case of a node experiment
    SASomaDist = SomaDist;
end

FPreStim = cell(NrExp,1);
Fbl      = cell(NrExp,1);
FdCal    = cell(NrExp,1);
FdFAP    = cell(NrExp,1);
FVm      = cell(NrExp,1);
FVmAbs   = cell(NrExp,1);
Decay    = exp(-SASomaDist/1180);

for iExp  = 1:NrExp
    FPreStim{iExp,1} = mean(FBCData{iExp}(:,PreStimFrames{iExp}),2,'omitnan');
    Fbl{iExp,1} = FBCData{iExp} - FPreStim{iExp};
    FdCal{iExp,1} = FCalInt{iExp} - FPreStim{iExp};
    DecayM = repmat(Decay,1,size(FBCData{iExp},2));
    FVm{iExp,1} = (Fbl{iExp} ./ FdCal{iExp}) .* dVm{iExp} .* DecayM;
    FVmAbs{iExp,1} = FVm{iExp} + VmRMP{iExp};
end

disp('Data is calibrated')

%% Filter data 
% I wrote a short script that uses a binomial filter with width of 3. 
% Define the data you want filtered (needs to be a cell) and the number of rounds

FFiltDataAbs = CFilter(FVmAbs,3);
FFiltData    = CFilter(FVm,3);

%% dVdt 
% Calculate the derivative of the membrane potential in the soma and axon.

%soma
VmdVdt = cellfun(@(x,q) diff(x)./diff(q), Vm, VmTime, 'UniformOutput', false);
[VmdVdtPeak, VmdVdtPeakInd] = cellfun(@(x) max(x,[],2), VmdVdt);
[VmdVdtMin, VmdVdtMinInd] = cellfun(@(x) min(x,[],2), VmdVdt);

for iExp = 1:NrExp   
    VmAPOnset(iExp,1) = VmTime{iExp}(VmdVdtPeakInd(iExp));
end

%axon
% To get AP Onset, I filter it 10 times
FFilt10xDataAbs = CFilter(FVmAbs,10);
FdVdt = cellfun(@(x,q) diff(x,[],2)./diff(q,[],2), FFilt10xDataAbs, FTime, 'UniformOutput',false);
[FdVdtPeak, FdVdtPeakInd] = cellfun(@(x) max(x,[],2), FdVdt, 'UniformOutput', false);
[FdVdtMin, FdVdtMinInd]   = cellfun(@(x) max(x,[],2), FdVdt, 'UniformOutput', false);

for iExp = 1:NrExp   
    for iPix = 1:29
        FAPOnset{iExp,1}(iPix,1) = FTime{iExp}(FdVdtPeakInd{iExp}(iPix));
    end
end

disp('dVdt is done') ;

%% Extract APPeak, APPeakInd, FStimStart, FTrialStart, SNR
% APPeak and APPeakInd are the peak and index of the peak of the AP.
% APPeak is taken from absolute Vm and FVm. APAmpTotal is relative from RMP
% VmStartStimMS is the start of the stimulus in ms for VmTime, is defined 
% DaMatlab.m. VmStimStart takes the corresponding acquisition point. It is
% only used to calculate the noise level before the stimulus starts.
% The Noise in VM and F is defined as the root mean square of the data
% acquired during baseline. SNR = (RMP - AP Peak) / noise

%soma
[VmAPPeak, VmAPPeakInd] = cellfun(@max, Vm);
VmAPAmpTotal            = VmAPPeak - cell2mat(VmRMP);
VmStimStart             = cellfun(@(x,q) find(x<(q-2),1,'last'), VmTime, VmStartStimMS, 'UniformOutput', false);   % -2 because StartStimMS can be inprecise, make sure to stop before onset stim. 
VmNoise                 = cellfun(@(x,idx,q) rms(x(:,idx)-q), Vm, VmStimStart, VmRMP, 'UniformOutput', false);
VmNoise                 = cell2mat(VmNoise);
VmSNR                   = VmAPAmpTotal ./ VmNoise;

%axon (peak within -1 and 1 ms in FTime
[trash, FPeakIndMinOne]     = cellfun(@(x) min(abs(x +1)), FTime, 'UniformOutput', false);
[trash, FPeakIndPlusOne]    = cellfun(@(x) min(abs(x -1)), FTime,'UniformOutput', false);
FStimStart                  = cellfun(@(x) x(end), PreStimFrames, 'UniformOutput', false);
[FAPAmpTotal, FAPPeakInd]   = cellfun(@(x,idx1,idx2) max(x(:,idx1:idx2),[],2), FFiltData, FPeakIndMinOne, FPeakIndPlusOne, 'UniformOutput', false);
FAPPeakInd                  = cellfun(@(x,q) x+q-1, FAPPeakInd, FPeakIndMinOne, 'UniformOutput', false);
FAPPeak                     = cellfun(@(x,q) x+q, FAPAmpTotal, VmRMP, 'UniformOutput', false);FTrialStart = cellfun(@(x) x(1), PreStimFrames, 'UniformOutput', false);

FNoise      = cellfun(@(x,idx) rms(x(:,idx),2), FVm, PreStimFrames, 'UniformOutput', false);
FSNR        = cellfun(@(x,q) x./q, FAPAmpTotal, FNoise,'UniformOutput', false);    

disp('Extracted APPeak, FTrialStart, FStimStart, SNR')

%% Threshold and AP Amplitude
% Threshold is defined as the membrane potential where dVdt passes 50 V/s
% In F, it is found by fitting a curve to the filtered data and then
% searching for the same steepness. APAmp is the amplitude of the AP
% relative from threshold.

%soma
VmThInd      = cellfun(@(x) find(x>50,1,'first'), VmdVdt, 'UniformOutput', false);
VmTh         = cellfun(@(x,idx) x(:,idx), Vm, VmThInd, 'UniformOutput', false);
VmThTime     = cellfun(@(x,idx) x(:,idx), VmTime, VmThInd, 'UniformOutput', false);
VmAPAmp      = VmAPAmpTotal - (cell2mat(VmTh) - cell2mat(VmRMP));

%axon (threshold found in filtered data)
ft = fittype('a*exp(-b*x) + c*exp(-d*x)');  %a and c determine amplitude peak, b and d more negative is sharper peak
fo = fitoptions('Method', 'NonlinearLeastSquares',...   %'Start',[12, -20, 15, -0.2]); %,...
    'Lower',[0,  -Inf, 0,   -Inf],...
    'Upper',[Inf, 0,   Inf, 0]);

FitFVm     = cell(y,NrExp);
FittedFVm  = cell(NrExp,1);
dFittedFVm = cell(NrExp,1);
FThFitInd  = cell(NrExp,1);
FTh        = cell(NrExp,1);
FThTime    = cell(NrExp,1);

for iExp = 1:NrExp
    for iPix = 1:y
        if ~isnan(FFiltData{iExp}(iPix,1))
            FitTo = FFiltData{iExp}(iPix,:)';
            xfit = FTime{iExp}';
            
            FitTo(1:FStimStart{iExp}) = nan;
            FitTo(FAPPeakInd{iExp}(iPix)-1:end) = nan;
            idxValid = ~isnan(FitTo);
                    
            xfit(1:FStimStart{iExp}) = nan;
            xfit(FAPPeakInd{iExp}(iPix):end) = nan;

            try
                FitFVm{iExp, iPix} = fit(xfit(idxValid), FitTo(idxValid), ft, fo);
                FittedFVm{iExp}(iPix,:) = feval(FitFVm{iExp, iPix},xfit)';
                dFittedFVm{iExp}(iPix,:) = diff(FittedFVm{iExp}(iPix,:)) ./ diff(xfit)';
                try
                    FThFitInd{iExp}(iPix,1) = find(dFittedFVm{iExp}(iPix,:) >50,1, 'first');
                    FTh{iExp}(iPix,1) = FittedFVm{iExp}(iPix,FThFitInd{iExp}(iPix));
                    FThTime{iExp}(iPix,1) = FTime{iExp}(FThFitInd{iExp}(iPix));
                catch
                    [~, FThFitInd{iExp}(iPix,1)] = max(dFittedFVm{iExp}(iPix,:));
                    FTh{iExp}(iPix,1) = FittedFVm{iExp}(iPix,FThFitInd{iExp}(iPix));
                    FThTime{iExp}(iPix,1) = FTime{iExp}(FThFitInd{iExp}(iPix));
                end
            catch
                FThFitInd{iExp}(iPix,1) = nan;
                FTh{iExp}(iPix,1) = nan;
                FThTime{iExp}(iPix,1) = nan;
            end
        else
        FThFitInd{iExp}(iPix,1) = nan;
        FTh{iExp}(iPix,1) = nan;
        FThTime{iExp}(iPix,1) = nan;
        end
    end
    disp(strcat('Threshold Progress:  ' , string(iExp), '/', string(NrExp)', ' done'));
end

FThAbs = cellfun(@(x,q) x+q, FTh, VmRMP, 'UniformOutput', false);
FAPAmp = cellfun(@(x,q,w) x - (q - w), FAPAmpTotal, FThAbs, VmRMP,'UniformOutput', false);

disp('Threshold and AP Amplitude is calculated')

%% FWHM
% FWHM or AP half width is found by the function fwhm.m. 
% Data is offset by the threshold. 
% Optionally you can fit the repolarization phase and replace the trace
% with that fit, to reduce the effect of noise. We tested this approach and
% the variance of the FWHM in our case did not reduce. (see supplemental
% figure)

% soma
VmBl   = cellfun(@(x,q) x-q, Vm, VmTh, 'UniformOutput', false);
VmFWHM = cellfun(@(x,q) fwhm(x,q), VmTime, VmBl, 'UniformOutput', false);

% axon
FFWHM = cell(NrExp,1);
%FFWHMfit = cell(NrExp,1);          

for iExp = 1:NrExp
    for iPix = 1:y
        try
        clear xtemp ytemp fitRange fittedY yToFit
            try 
                xtemp = FTime{iExp}(FThFitInd{iExp}(iPix)-2:FAPPeakInd{iExp}(iPix)+100);   
                ytemp = FFiltData{iExp}(iPix,FThFitInd{iExp}(iPix)-2:FAPPeakInd{iExp}(iPix)+100);    
            catch
                xtemp = FTime{iExp}(FThFitInd{iExp}(iPix)-2:end);   
                ytemp = FFiltData{iExp}(iPix,FThFitInd{iExp}(iPix)-2:end);    
            end

        ytemp = ytemp - FTh{iExp}(iPix);
        
%         fitRange = (FAPPeakInd{iExp}(iPix)-FThFitInd{iExp}(iPix)+2):(FAPPeakInd{iExp}(iPix)-FThFitInd{iExp}(iPix)+40);
%         fitY = fit(xtemp(fitRange)', ytemp(fitRange)', 'Poly2');     
%         fittedY = feval(fitY, xtemp(fitRange));
%         yToFit = ytemp;
%         yToFit(fitRange) = fittedY;
       
        FFWHM{iExp}(iPix,1) = fwhm(xtemp, ytemp);
        %FFWHMfit{iExp}(iPix,1) = fwhm(xtemp, yToFit);
        
        catch
        FFWHM{iExp}(iPix,1) = NaN;
        %FWHMfit{iExp}(iPix,1) = NaN;
        end
    end
    disp(strcat('FWHM Progress:  ' , string(iExp), '/', string(NrExp)', ' done'));
end


disp('FWHM is calculated')

%% AHP and ADP
% In control patch-clamp, AHP is defined as the minumun after the AP peak 
% and ADP as the small maximum afer that minimum. 
% For non-ctrl conditions, instead take the value at the index of the AHP
% and ADP in ctrl. 
% Disadvantage is the rare case that ADP is badly measured in ctrl (for
% example when it was during pipette artifact), and then it is taken at the
% wrong time in other conditions. 
% In the axon, use the time of the Vm AHP and ADP to fit around that moment
% AHP and ADP are found by fitting curves in F and then I find the local
% dip and maxima in ctrl. I use the control index in non-ctrl recordings to
% find amplitude at that time, exactly like in the patch-clamp recording.
% First, a curve is fitted to the fluorescnece data from the AP Peak until
% 20 frames (=1ms) after the AHP index in the soma. A Poly4 fit is used.
% You need at least a 3rd polynomial to fit the steepness of the AP in 
% combination with the dip plus another rise.
% See supplemental figure X for justification of that fit. The dip in this 
% fit is defined as the AHP. Subsequently, a second curve is fitted from
% the AHP index until the end of the trace. The maximum of this trace is
% the ADP. 

% Soma
VmAHPRangeS = num2cell([VmAPPeakInd+20]);
VmAHPRangeE = num2cell([VmAPPeakInd+300]);

[VmAHP, VmAHPInd]  = cellfun(@(x,idxS,idxE) min(x(:,idxS:idxE)), VmBl, VmAHPRangeS, VmAHPRangeE, 'UniformOutput', false);
VmAHPInd           = cellfun(@(x,q) x+q-1, VmAHPInd, VmAHPRangeS, 'UniformOutput', false);

[VmADP, VmADPInd]  = cellfun(@(x,idx) max(x(:,idx:end)), VmBl, VmAHPInd, 'UniformOutput', false);
VmADPInd           = cellfun(@(x,q) x+q-1, VmADPInd, VmAHPInd, 'UniformOutput', false);

cellIdCat          = categorical(cellId);
conditionCat       = categorical(condition);

% for iExp = 1:NrExp
%     if conditionCat(iExp) == 'ctrl'         
%     else
%         tempCell = cellId{iExp};
%         [val,~] = intersect(find(strcmp(condition,'ctrl')==1),find(cellIdCat == tempCell));
%         VmAHPInd{iExp} = VmThInd{iExp}+(VmAHPInd{val} - VmThInd{val});
%         VmADPInd{iExp} = VmThInd{iExp}+(VmADPInd{val} - VmThInd{val});
%         VmAHP{iExp} = VmBl{iExp}(VmAHPInd{iExp});
%         VmADP{iExp} = VmBl{iExp}(VmADPInd{iExp});
%     end
% end

VmAHPTime = cellfun(@(x,y) x(y), VmTime, VmAHPInd, 'UniformOutput', false);
VmADPTime = cellfun(@(x,y) x(y), VmTime, VmADPInd, 'UniformOutput', false);

% Axon.
xfit        = cell(NrExp, y);
FAHPFit     = cell(NrExp, y);
FAHPFitted  = cell(NrExp, y);
FAHPrmp     = cell(NrExp, 1);
FAHPTime    = cell(NrExp, 1);
FAHPInd     = cell(NrExp, 1);
FAHPxfit    = cell(NrExp, 1);

FADPFit     = cell(NrExp, y);
FADPFitted  = cell(NrExp, y);
FADPrmp     = cell(NrExp, 1);
FADPTime    = cell(NrExp, 1);           
FADPInd     = cell(NrExp, 1);
FADPxfit    = cell(NrExp, 1);

for iExp = 1:NrExp
    if conditionCat(iExp) == 'ctrl'                                                 % Find AHP and ADP in ctrl first. 
        for iPix = 1:y
            [~,closestIndexAHP] = min(abs(VmAHPTime{iExp}-FTime{iExp}));            % find in FTime closest value to AHP index 
                 
            xfit  = FTime{iExp}(1,FAPPeakInd{iExp}(iPix)+1:closestIndexAHP+20);             
            FitTo = FFiltData{iExp}(iPix,FAPPeakInd{iExp}(iPix)+1:closestIndexAHP+20);    
                
            if ~isnan(FitTo)
                try
                    
                    FAHPFit{iExp,iPix}    = fit(xfit',FitTo','Poly4');     
                    FAHPFitted{iExp,iPix} = feval(FAHPFit{iExp,iPix} ,xfit);
                                        
                    tempLog = islocalmin(FAHPFitted{iExp,iPix});
                    tempAHP = FAHPFitted{iExp,iPix}(tempLog);
                    tempAHPTime = xfit(tempLog);
                    tempAHP(2:end) = [];
                    tempAHPTime(2:end) = [];
                                        
                    if isempty(tempAHP)
                        [tempAHP, tempAHPInd] = min(FAHPFitted{iExp,iPix});
                        tempAHPTime = xfit(tempAHPInd);
                    end                     
                    FAHPrmp{iExp}(iPix,1) = tempAHP;
                    FAHPTime{iExp}(iPix,1) = tempAHPTime;
                    FAHPInd{iExp}(iPix,1) = find(FTime{iExp} == tempAHPTime); 
                    FAHPxfit{iExp,iPix} = xfit; 

                    xfit  = FTime{iExp}(1,FAHPInd{iExp}(iPix):end);
                    FitTo = FFiltData{iExp}(iPix,FAHPInd{iExp}(iPix):end); 

                    FADPFit{iExp,iPix}    = fit(xfit',FitTo','Poly4');     
                    FADPFitted{iExp,iPix} = feval(FADPFit{iExp,iPix} ,xfit);
                    
                    tempLog = islocalmax(FADPFitted{iExp,iPix});
                    tempADP = FADPFitted{iExp,iPix}(tempLog);
                    tempADPTime = xfit(tempLog);
                    tempADP(2:end) = [];
                    tempADPTime(2:end) = [];
                    
                    if isempty(tempADP)
                        [tempADP, tempADPInd] = max(FADPFitted{iExp,iPix});
                        tempADPTime = xfit(tempADPInd);
                    end 
                    FADPrmp{iExp}(iPix,1) = tempADP;
                    FADPTime{iExp}(iPix,1) = tempADPTime;
                    FADPInd{iExp}(iPix,1) = find(FTime{iExp} == tempADPTime); 
                    FADPxfit{iExp,iPix} = xfit;
                    clear FitTo tempLog tempAHP tempAHPTime tempADP tempADPTime
                catch
                    FAHPrmp{iExp}(iPix,1) = nan;
                    FAHPTime{iExp}(iPix,1) = nan;
                    FADPrmp{iExp}(iPix,1) = nan;
                    FADPTime{iExp}(iPix,1) = nan;
                end
            else
                FAHPrmp{iExp}(iPix,1) = nan;
                FAHPTime{iExp}(iPix,1) = nan;
                FADPrmp{iExp}(iPix,1) = nan;
                FADPTime{iExp}(iPix,1) = nan;
            end
        end
    end
    disp(strcat('AHP Progress (ctrl):  ' , string(iExp), '/', string(NrExp)', ' done'));
end

% AHP and ADP in control are found. Now use that index to find them in
% non-control condition. 

for iExp = 1:NrExp
    if conditionCat(iExp) == 'ctrl'
    else
        tempCell = cellId{iExp};
        for iPix = 1:y
            [val,~] = intersect(find(strcmp(condition,'ctrl')==1),find(cellIdCat == tempCell));
            [~, FAHPInd{iExp}(iPix,1)] = min(abs(FTime{iExp} - FAHPTime{val}(iPix,1))); %absolute time doesnt always match, so find min difference
            [~, FADPInd{iExp}(iPix,1)] = min(abs(FTime{iExp} - FADPTime{val}(iPix,1)));   
            try         % try because in some rare cases the index + 20 exceeds the data
                xfit  = FTime{iExp}(1,FAPPeakInd{iExp}(iPix)+1:FAHPInd{iExp}(iPix,1)+20); 
                FitTo = FFiltData{iExp}(iPix,FAPPeakInd{iExp}(iPix)+1:FAHPInd{iExp}(iPix,1)+20); 
            catch        
                xfit  = FTime{iExp}(1,FAPPeakInd{iExp}(iPix)+1:end); 
                FitTo = FFiltData{iExp}(iPix,FAPPeakInd{iExp}(iPix)+1:end);
            end
        
            if ~isnan(FitTo)
                try
                    FAHPFit{iExp,iPix}    = fit(xfit',FitTo','Poly4');         % tried many types of fits, this seems to work best. 
                    FAHPFitted{iExp,iPix} = feval(FAHPFit{iExp,iPix} ,xfit);
                                        
                    FAHPrmp{iExp}(iPix,1) = FAHPFitted{iExp,iPix}((FAHPInd{iExp}(iPix)+1 - FAPPeakInd{iExp}(iPix)+1),1);
                    FAHPTime{iExp}(iPix,1) = FTime{iExp}(FAHPInd{iExp}(iPix,1));
                    FAHPxfit{iExp,iPix} = xfit; 
                catch
                    FAHPrmp{iExp}(iPix,1) = nan;
                    FAHPTime{iExp}(iPix,1) = nan;
                end
                
                try
                    xfit  = FTime{iExp}(1,FAHPInd{iExp}(iPix):end);
                    FitTo = FFiltData{iExp}(iPix,FAHPInd{iExp}(iPix):end); 
                   
                    FADPFit{iExp,iPix}    = fit(xfit',FitTo','Poly4');     % tried many types, this seems to work best.
                    FADPFitted{iExp,iPix} = feval(FADPFit{iExp,iPix} ,xfit);
                     
                    FADPrmp{iExp}(iPix,1) = FADPFitted{iExp,iPix}((FADPInd{iExp}(iPix)+1 - FAHPInd{iExp}(iPix)),1);  
                    FADPTime{iExp}(iPix,1) = FTime{iExp}(FADPInd{iExp}(iPix,1));
                    FADPxfit{iExp,iPix} = xfit;                        
                catch
                    FADPrmp{iExp}(iPix,1) = nan;
                    FADPTime{iExp}(iPix,1) = nan; 
                end
                clear FitTo tempLog tempAHP tempAHPTime tempADP tempADPTime                 
      
            else
                FAHPrmp{iExp}(iPix,1) = nan;
                FAHPTime{iExp}(iPix,1) = nan;
                FADPrmp{iExp}(iPix,1) = nan;
                FADPTime{iExp}(iPix,1) = nan;
            end
        end
    end
    disp(strcat('AHP Progress (non-ctrl):  ' , string(iExp), '/', string(NrExp)', ' done'));
end

FAHP = cellfun(@(x,q) x-q, FAHPrmp, FTh,'UniformOutput', false);
FADP = cellfun(@(x,q) x-q, FADPrmp, FTh,'UniformOutput', false);

disp('AHP and ADP are calculated')

%% Area under the curve
% I also quantify the area under the curve as a measure of repolarization
% that is less affected by SNR, but more by calibration.

% Soma

VmZero = cellfun(@(x,q) x-q, Vm, VmRMP, 'UniformOutput', false);

for iExp = 1:NrExp
    
    VmAUCy    = VmZero{iExp}(VmThInd{iExp}:VmAHPInd{iExp});
    AUCtemp = nan(size(VmTime{iExp}));
    AUCtemp((VmThInd{iExp}:VmAHPInd{iExp})) = cumtrapz(VmAUCy)/100;
        
    VmAUCcum{iExp,1} = AUCtemp;
    VmAUC{iExp,1}    = max(AUCtemp);
    
end

%axon

FVmZero = cellfun(@(x,q) x-q, FVmAbs, VmRMP, 'UniformOutput', false);
FAUCcum = cell(NrExp, 1);
FAUC    = cell(NrExp, 1);

for iExp = 1:NrExp
    for iPix = 1:y
        try
            ytemp = FVmZero{iExp}(iPix,FThFitInd{iExp}(iPix):FAHPInd{iExp}(iPix));
            AUCtemp = nan(size(FTime{iExp}));
            AUCtemp(FThFitInd{iExp}(iPix):FAHPInd{iExp}(iPix)) = cumtrapz(ytemp)/20;
        
            FAUCcum{iExp}(iPix,:) = AUCtemp;
            FAUC{iExp}(iPix,1) = max(AUCtemp);
            
        catch
            FAUCcum{iExp}(iPix,:) = nan;
            FAUC{iExp}(iPix,1) = nan;
        end
    end
end

disp('AUC is calculated')

%% Now combine Ephys and imaging parameters

APAmp   = [VmAPAmp             cell2mat(cellfun(@(x) transpose(x), FAPAmp, 'UniformOutput', false))]; 
APOnset = [VmAPOnset           cell2mat(cellfun(@(x) transpose(x), FAPOnset, 'UniformOutput', false))];    
dVdtP   = [VmdVdtPeak          cell2mat(cellfun(@(x) transpose(x), FdVdtPeak, 'UniformOutput', false))];
dVdtPInd= [VmdVdtPeakInd       cell2mat(cellfun(@(x) transpose(x), FdVdtPeakInd, 'UniformOutput', false))];
dVdtM   = [VmdVdtMin           cell2mat(cellfun(@(x) transpose(x), FdVdtMin, 'UniformOutput', false))]; 
Th      = [cell2mat(VmTh)      cell2mat(cellfun(@(x) transpose(x), FThAbs, 'UniformOutput', false));];
ThTime  = [cell2mat(VmThTime)  cell2mat(cellfun(@(x) transpose(x), FThTime, 'UniformOutput', false))];
FWHM    = [cell2mat(VmFWHM)    cell2mat(cellfun(@(x) transpose(x), FFWHM, 'UniformOutput', false))];
%FWHMfit = [cell2mat(VmFWHM)   cell2mat(cellfun(@(x) transpose(x), FFWHMfit, 'UniformOutput', false))];
AHP     = [cell2mat(VmAHP)     cell2mat(cellfun(@(x) transpose(x), FAHP, 'UniformOutput', false))];
ADP     = [cell2mat(VmADP)     cell2mat(cellfun(@(x) transpose(x), FADP, 'UniformOutput', false))];
%AUC    = [cell2mat(VmAUC)     cell2mat(cellfun(@(x) transpose(x), FAUC, 'UniformOutput', false))];
SNR     = [VmSNR               cell2mat(cellfun(@(x) transpose(x), FSNR, 'UniformOutput', false))]; 
AHPTime = [cell2mat(VmAHPTime) cell2mat(cellfun(@(x) transpose(x), FAHPTime, 'UniformOutput', false))]; 


%% Discard pixels with SNR < 2 and very obvious outliers

APAmp(SNR <2) = nan;
Th(SNR <2) = nan;
FWHM(SNR <2) = nan;
ThTime(SNR <2) = nan;
AHP(SNR <2) = nan;
AHPTime(SNR <2) = nan;

% exclusion criteria
ThTime(ThTime <-0.5) = nan;
ThTime(ThTime >1) = nan;

FWHMNanBefore = sum(sum(isnan(FWHM)));
FWHM(FWHM>1.75) = NaN;            % Widest AP is after 4-AP application is 1.18 ms, so bigger than 1.75 must be outlier
FWHM(FWHM == 0) = NaN;
FWHMDiscarded = sum(sum(isnan(FWHM))) - FWHMNanBefore;


%% Control data
% calculate mean and sem for given parameter to plot

condition = string(condition);

Param = Th(condition == 'ctrl',:);
ParamMean = mean(Param,1, 'omitnan');
ParaSD = std(Param, 'omitnan');
ParamNumber = sum(~isnan(Param),1);
ParaSem = ParaSD ./ sqrt(ParamNumber);
ParamMeanSem = [(ParamMean - ParaSem)' ParamMean' (ParamMean + ParaSem)'];
ParaSomaAisInt = [Param(:,1) Param(:,28:30)];
ParamMedian = median(Param,1, 'omitnan');
ParamMedianSem = [(ParamMedian - ParaSem)' ParamMedian' (ParamMedian + ParaSem)'];


%% pair per condition
% Adjust condition and parameter (ParToPlot) to get the data you can copy
% to Prism for a grouped comparison

condition = string(condition);
clear DataPrism PrismC1 PrismC2 PrismC3
ParToPlot = dVdtP;
ConToPlot = {'ctrl','cavs','cavs'};
CellsToPlot = cellId(condition == ConToPlot{3});
[~, locb] = ismember(cellId, CellsToPlot);

% Index conditions
cond1 = ConToPlot{1};
cond2 = ConToPlot{2};
cond3 = ConToPlot{3};
NrCond = length(ConToPlot);

for i = 1:length(cellId)
    if condition{i} == cond1
        cond1Ind(i) = 1;
    elseif condition{i} == cond2
        cond2Ind(i) = 1;
    elseif condition{i} == cond3
        cond3Ind(i) = 1; 
    end
end

% Index cells
Cells = cellId(1);
cellId = cellstr(cellId);
for teller = 2:length(cellId)
    if cellId{teller} == cellId{teller-1}
    else
        Cells(teller) = cellId(teller);
    end
end
Cells = Cells(~cellfun('isempty',Cells));


cellIdCat = categorical(cellId);
NrCellsToPair = length(CellsToPlot);

% Now create empty matrix that can be filled with data in the correct
% order. Also leaves blanks for conditions that weren't recorded.

if NrCond == 2
    CondOrder = repmat([{cond1}; {cond2}],NrCellsToPair,1);
elseif NrCond == 3
    CondOrder = repmat([{cond1}; {cond2}; {cond3}],NrCellsToPair,1);
end

% first create matrix with conditions in the right order
Data = nan(NrCellsToPair*NrCond,y+1);
DataN = nan(NrCellsToPair*NrCond,y+1);
% Also create empty matrix that can be filled with traces from the selected
% pixel as examples
DataToPlot = nan(NrCellsToPair*NrCond,y+1, 4);


for iCell = 1:NrCellsToPair
    iCond = 1;
    cellToWrite = CellsToPlot(iCell);
    condToWrite = CondOrder(iCond); 
    [val,~] = intersect(find(strcmp(condition,condToWrite)==1),find(cellIdCat == cellToWrite));
    if(~isnan(val))
        A  = ParToPlot(val,:);
    else
        A = nan(1,y+1);
    end
    iCond = iCond + 1;
    condToWrite = CondOrder(iCond);   
    [val,~] = intersect(find(strcmp(condition,condToWrite)==1),find(cellIdCat == cellToWrite));
    if(~isnan(val))
        B  = ParToPlot(val,:);        
    else
        B = nan(1,y+1);
    end
    PrismC1(iCell,:) = A;
    PrismC2(iCell,:) = B;
    if NrCond == 2
        Data(((NrCond*iCell)- NrCond +1):((NrCond*iCell)- NrCond +2),:) = [A; B];
        DataN(((NrCond*iCell)- NrCond +1):((NrCond*iCell)- NrCond +2),:) = [A./A; B./A];
        
        
    elseif NrCond == 3
        iCond = iCond + 1;
        condToWrite = CondOrder(iCond);   
        [val,~] = intersect(find(strcmp(condition,condToWrite)==1),find(cellIdCat == cellToWrite));
        if(~isnan(val))
            C  = ParToPlot(val,:);         
        else
            C = nan(1,y+1);        
        end
        Data(((NrCond*iCell)- NrCond +1):((NrCond*iCell)- NrCond +3),:) = [A; B; C];
        DataN(((NrCond*iCell)- NrCond +1):((NrCond*iCell)- NrCond +3),:) = [A./A; B./A; C./A];
        PrismC3(iCell,:) = C;
    end
    clear A B C
    %disp(iCell)
end


if NrCond == 3
    DataPrism = [PrismC1(:,1)' PrismC2(:,1)' PrismC3(:,1)'; PrismC1(:,28:30)'  PrismC2(:,28:30)' PrismC3(:,28:30)'];
elseif NrCond == 2
    DataPrism = [PrismC1(:,1)' PrismC2(:,1)'; PrismC1(:,28:30)'  PrismC2(:,28:30)'];
end


%% Plot traces from Soma, Proximal AIS, distal AIS and internode for each individual recording
% Including all fits and estimated parameters



for iExp = 1:NrExp
    figure('Position', [800 10 500 1300])
    cellToPlot = cellId{iExp};
    PlottedCond = condition{iExp}; 
    subplot(4,1,1)
    hold on
    title(strcat(cellToPlot(1:6),'-', cellToPlot(8:9),'-', PlottedCond))
    plot(VmTime{iExp}, Vm{iExp})
    plot(VmTime{iExp}(VmThInd{iExp}),Vm{iExp}(VmThInd{iExp}),'rx')
    plot(VmTime{iExp}(VmAHPInd{iExp}),VmAHP{iExp}+VmTh{iExp},'rx')  
    plot(VmTime{iExp}(VmADPInd{iExp}),VmADP{iExp}+VmTh{iExp},'ro')
    plot([7 7],[-25 25], 'k-')
    plot([6 7],[-25 -25], 'k-')
    %xlim([-10 ])
    ylim([-100 70])
 
    subplot(4,1,2) 
    hold on
    plot(FTime{iExp}, FFiltDataAbs{iExp}(27,:))
    plot(FTime{iExp}, FittedFVm{iExp}(27,:)+VmRMP{iExp},'k')
    plot(FTime{iExp}(FThFitInd{iExp}(27,:)), FTh{iExp}(27,:)+VmRMP{iExp}, 'rx') 
    try       % not always an AHP, so try
        plot(FAHPxfit{iExp,27}, FAHPFitted{iExp,27}(:)+VmRMP{iExp})
        plot(FAHPTime{iExp}(27), FAHPrmp{iExp}(27)+VmRMP{iExp},'rx')
    catch
    end
    try     % not always an ADP, so try
        plot(FADPxfit{iExp,27}, FADPFitted{iExp,27}(:)+VmRMP{iExp})
        plot(FADPTime{iExp}(27), FADPrmp{iExp}(27)+VmRMP{iExp},'ro') 
    catch
    end
    title(strcat('SNR: ', string(SNR(iExp,28))))         % One higher because SNR matrix contains Soma data in first column
    %xlim([-4 6])
    ylim([-100 70])
   
    subplot(4,1,3)
    hold on
    plot(FTime{iExp}, FFiltDataAbs{iExp}(28,:))
    plot(FTime{iExp}, FittedFVm{iExp}(28,:)+VmRMP{iExp},'k')
    plot(FTime{iExp}(FThFitInd{iExp}(28,:)), FTh{iExp}(28,:)+VmRMP{iExp}, 'rx')
    try
        plot(FAHPxfit{iExp,28}, FAHPFitted{iExp,28}(:)+VmRMP{iExp})
        plot(FAHPTime{iExp}(28), FAHPrmp{iExp}(28)+VmRMP{iExp},'rx')
    catch
    end
    try
        plot(FADPxfit{iExp,28}, FADPFitted{iExp,28}(:)+VmRMP{iExp})
        plot(FADPTime{iExp}(28), FADPrmp{iExp}(28)+VmRMP{iExp},'ro')    
    catch
    end 
    title(strcat('SNR: ', string(SNR(iExp,29))))
    %xlim([-4 6])
    ylim([-100 70])

    subplot(4,1,4)
    hold on
    plot(FTime{iExp}, FFiltDataAbs{iExp}(29,:))
    plot(FTime{iExp}, FittedFVm{iExp}(29,:)+VmRMP{iExp},'k')
    plot(FTime{iExp}(FThFitInd{iExp}(29,:)), FTh{iExp}(29,:)+VmRMP{iExp}, 'rx')
    try
        plot(FAHPxfit{iExp,29}, FAHPFitted{iExp,29}(:)+VmRMP{iExp})
        plot(FAHPTime{iExp}(29), FAHPrmp{iExp}(29)+VmRMP{iExp},'rx')
    catch
    end
    try
        plot(FADPxfit{iExp,29}, FADPFitted{iExp,29}(:)+VmRMP{iExp})  
        plot(FADPTime{iExp}(29), FADPrmp{iExp}(29)+VmRMP{iExp},'ro')   
    catch
    end   
    title(strcat('SNR: ', string(SNR(iExp,30))))
    %xlim([-4 6])
    ylim([-100 70])
 
    saveas(gcf,strcat(outputpath, 'Figures/MatlabExports/SingleCellConRecCGH/',string(cellId{iExp}), '_' ,string(condition{iExp})))
    close
    disp(strcat('Progress:  ' , string(iExp), '/', string(NrExp)', ' done'));
end


%% plot example traces from paired data
% aligned on RMP

UniqueCells = unique(cellId);

% [1 3 2];

for iCell = 2:length(UniqueCells)
    figure('Position', [800 10 500 1300])
    cellToPlot = UniqueCells{iCell};
    
    ind = (cellId == string(cellToPlot));
    condsToPlot = condition(ind);     
    ExpsToPlot = (1:1:NrExp);
    ExpsToPlot = ExpsToPlot(ind);
    NrCond = size(condsToPlot,1);
    
    if NrCond == 1 
        colorOrder = {'k'};
        savename = strcat(string(cellToPlot), '_' ,string(condsToPlot(1)), '_PairedCondTraces');
    elseif NrCond == 2
        if condsToPlot(1) == "hg09" 
        else
        ExpsToPlot = ExpsToPlot([2 1]);
        condsToPlot = condsToPlot([2 1]) ;
        end
        colorOrder = {'k', 'r'};   
        savename = strcat(string(cellToPlot), '_' ,string(condsToPlot(2)), '_PairedCondTraces');
    elseif NrCond == 3        
        if     condsToPlot == ["20mi"; "40mi"; "ctrl"]        
            CondOrder = [3 1 2];
        elseif condsToPlot == ["bapt"; "bpIb"; "ctrl"]  
            CondOrder = [3 1 2];
        elseif condsToPlot == ["ctrl"; "egIb"; "egta"]
            CondOrder = [1 3 2]; 
        elseif condsToPlot == ["ctrl"; "px01"; "px10"]  
            CondOrder = [1 2 3];
        elseif condsToPlot == ["4-AP"; "4Ibx"; "ctrl"]
            CondOrder = [3 1 2];
        elseif condsToPlot == ["cavs"; "ctrl"; "cvIb"]  
            CondOrder = [1 2 3];
        elseif condsToPlot == ["ctrl"; "ryIb"; "ryan"]
            CondOrder = [1 3 2];
        elseif condsToPlot == ["hg09"; "hg10"; "squa"]
            CondOrder = [3 2 1];
        else  
            disp(condsToPlot)
            CondOrder = input('What is the right order?') 
        end

        ExpsToPlot = ExpsToPlot(CondOrder);
        condsToPlot = condsToPlot(CondOrder);
        colorOrder = {'k', 'r', 'b'};
        savename = strcat(string(cellToPlot), '_' ,string(condsToPlot(2)),'_' ,string(condsToPlot(3)), '_PairedCondTraces');
    end
  
    for iCond = 1:NrCond
        subplot(4,1,1)
        title(strcat(cellToPlot(1:6),'-', cellToPlot(8:9)))
        hold on
        plot(VmTime{ExpsToPlot(iCond)}, Vm{ExpsToPlot(iCond)},colorOrder{iCond})
        %plot(VmTime{ExpsToPlot(iCond)}(VmAHPInd{ExpsToPlot(iCond)}),VmAHP{ExpsToPlot(iCond)}+VmTh{ExpsToPlot(iCond)},'rx')  
        %plot(VmTime{ExpsToPlot(iCond)}(VmADPInd{ExpsToPlot(iCond)}),VmADP{ExpsToPlot(iCond)}+VmTh{ExpsToPlot(iCond)},'ro')
        plot([5.5 5.5],[-25 25], 'k-')
        plot([4.5 5.5],[-25 -25], 'k-')
        xlim([-3 6])
        ylim([-100 70])
        legend(condsToPlot); 
    end
    
    for iCond = 1:NrCond
        subplot(4,1,2)
        hold on
        plot(FTime{ExpsToPlot(iCond)}, FFiltDataAbs{ExpsToPlot(iCond)}(27,:),colorOrder{iCond})
%         plot(FAHPxfit{ExpsToPlot(iCond),27}, FAHPFitted{ExpsToPlot(iCond),27}(:)+VmRMP{ExpsToPlot(iCond)}, colorOrder{iCond})
%         plot(FAHPTime{ExpsToPlot(iCond)}(27), FAHPrmp{ExpsToPlot(iCond)}(27)+VmRMP{ExpsToPlot(iCond)},'rx')
%         plot(FADPxfit{ExpsToPlot(iCond),27}, FADPFitted{ExpsToPlot(iCond),27}(:)+VmRMP{ExpsToPlot(iCond)}, colorOrder{iCond})
%         plot(FADPTime{ExpsToPlot(iCond)}(27), FADPrmp{ExpsToPlot(iCond)}(27)+VmRMP{ExpsToPlot(iCond)},'ro')
        xlim([-3 6])
        ylim([-100 70])
        LegendToBe(iCond,:) = strcat('SNR:',string(SNR(ExpsToPlot(iCond),28)));  
    end
    legend(LegendToBe)
    
    for iCond = 1:NrCond
        subplot(4,1,3)
        hold on
        plot(FTime{ExpsToPlot(iCond)}, FFiltDataAbs{ExpsToPlot(iCond)}(28,:),colorOrder{iCond})
%         plot(FAHPxfit{ExpsToPlot(iCond),28}, FAHPFitted{ExpsToPlot(iCond),28}(:)+VmRMP{ExpsToPlot(iCond)}, colorOrder{iCond})
%         plot(FAHPTime{ExpsToPlot(iCond)}(28), FAHPrmp{ExpsToPlot(iCond)}(28)+VmRMP{ExpsToPlot(iCond)},'rx')
%         plot(FADPxfit{ExpsToPlot(iCond),28}, FADPFitted{ExpsToPlot(iCond),28}(:)+VmRMP{ExpsToPlot(iCond)}, colorOrder{iCond})
%         plot(FADPTime{ExpsToPlot(iCond)}(28), FADPrmp{ExpsToPlot(iCond)}(28)+VmRMP{ExpsToPlot(iCond)},'ro')
        xlim([-3 6])
        ylim([-100 70])
        LegendToBe(iCond,:) = strcat('SNR:',string(SNR(ExpsToPlot(iCond),29)));  
    end
    legend(LegendToBe)
 
    for iCond = 1:NrCond
        subplot(4,1,4)
        hold on
        plot(FTime{ExpsToPlot(iCond)}, FFiltDataAbs{ExpsToPlot(iCond)}(29,:),colorOrder{iCond})
%         plot(FAHPxfit{ExpsToPlot(iCond),29}, FAHPFitted{ExpsToPlot(iCond),29}(:)+VmRMP{ExpsToPlot(iCond)}, colorOrder{iCond})
%         plot(FAHPTime{ExpsToPlot(iCond)}(29), FAHPrmp{ExpsToPlot(iCond)}(29)+VmRMP{ExpsToPlot(iCond)},'rx')
%         plot(FADPxfit{ExpsToPlot(iCond),29}, FADPFitted{ExpsToPlot(iCond),29}(:)+VmRMP{ExpsToPlot(iCond)}, colorOrder{iCond})
%         plot(FADPTime{ExpsToPlot(iCond)}(29), FADPrmp{ExpsToPlot(iCond)}(29)+VmRMP{ExpsToPlot(iCond)},'ro')
        xlim([-3 6])
        ylim([-100 70])
        LegendToBe(iCond,:) = strcat('SNR:',string(SNR(ExpsToPlot(iCond),30)));  
    end
    legend(LegendToBe)
 
    saveas(gcf,strcat(outputpath, 'Figures/MatlabExports/PairedTraces/',savename))
    %pause
    close
    disp(strcat('Progress:  ' , string(iCell), '/', string(length(UniqueCells))', ' done'));
end
%% plot all traces from aligned paired data
% aligned on somatic threshold

UniqueCells = unique(cellId);


for iCell = 1:length(UniqueCells)
    figure('Position', [800 10 500 1300])
    cellToPlot = UniqueCells{iCell};
    
    ind = (cellId == string(cellToPlot));
    condsToPlot = condition(ind);     
    ExpsToPlot = (1:1:NrExp);
    ExpsToPlot = ExpsToPlot(ind);
    NrCond = size(condsToPlot,1);
    
    if NrCond == 1 
        colorOrder = {'k'};
        savename = strcat(string(cellToPlot), '_' ,string(condsToPlot(1)), '_PairedAlignedCondTraces');
    elseif NrCond == 2
        if condsToPlot(1) == 'ctrl'
        else
        ExpsToPlot = ExpsToPlot([2 1]);
        condsToPlot = condsToPlot([2 1]) ;
        end
        colorOrder = {'k', 'r'};   
        savename = strcat(string(cellToPlot), '_' ,string(condsToPlot(2)), '_PairedAlignedCondTraces');
    elseif NrCond == 3        
        if     condsToPlot == ["20mi"; "40mi"; "ctrl"]        
            CondOrder = [3 1 2];
        elseif condsToPlot == ["bapt"; "bpIb"; "ctrl"]  
            CondOrder = [3 1 2];
        elseif condsToPlot == ["ctrl"; "egIb"; "egta"]
            CondOrder = [1 3 2]; 
        elseif condsToPlot == ["ctrl"; "px01"; "px10"]  
            CondOrder = [1 2 3];
        elseif condsToPlot == ["4-AP"; "4Ibx"; "ctrl"]
            CondOrder = [3 1 2];
        elseif condsToPlot == ["cavs"; "ctrl"; "cvIb"]  
            CondOrder = [1 2 3];
        elseif condsToPlot == ["ctrl"; "ryIb"; "ryan"]
            CondOrder = [1 3 2];
        elseif condsToPlot == ["hg09"; "hg10"; "squa"]
            CondOrder = [3 2 1];
        else  
            disp(condsToPlot)
            CondOrder = input('What is the right order?') 
        end
   
        ExpsToPlot = ExpsToPlot(CondOrder);
        condsToPlot = condsToPlot(CondOrder);
        colorOrder = {'k', 'r', 'b'};
        savename = strcat(string(cellToPlot), '_' ,string(condsToPlot(2)),'_' ,string(condsToPlot(3)), '_PairedAlignedCondTraces');
    end
  
    for iCond = 1:NrCond
        subplot(4,1,1)
        title(strcat(cellToPlot(1:6),'-', cellToPlot(8:9)))
        hold on
        plot(VmTime{ExpsToPlot(iCond)}, Vm{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)},colorOrder{iCond})
        plot(VmTime{ExpsToPlot(iCond)}(VmAHPInd{ExpsToPlot(iCond)}),VmAHP{ExpsToPlot(iCond)},'rx')  
        plot(VmTime{ExpsToPlot(iCond)}(VmADPInd{ExpsToPlot(iCond)}),VmADP{ExpsToPlot(iCond)},'ro')    
        plot([5.5 5.5],[25 75], 'k-')
        plot([4.5 5.5],[25 25], 'k-')
        xlim([-3 6])
        ylim([-50 140])
        legend(condsToPlot); 
    end
    
    for iCond = 1:NrCond
        subplot(4,1,2)
        hold on
        plot(FTime{ExpsToPlot(iCond)}, FFiltDataAbs{ExpsToPlot(iCond)}(27,:)-VmTh{ExpsToPlot(iCond)},colorOrder{iCond})
        plot(FAHPxfit{ExpsToPlot(iCond),27}, FAHPFitted{ExpsToPlot(iCond),27}(:)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)}, colorOrder{iCond})
        plot(FAHPTime{ExpsToPlot(iCond)}(27), FAHPrmp{ExpsToPlot(iCond)}(27)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)},'rx')
        plot(FADPxfit{ExpsToPlot(iCond),27}, FADPFitted{ExpsToPlot(iCond),27}(:)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)}, colorOrder{iCond})
        plot(FADPTime{ExpsToPlot(iCond)}(27), FADPrmp{ExpsToPlot(iCond)}(27)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)},'ro')
        xlim([-3 6])
        ylim([-50 140])
        LegendToBe(iCond,:) = strcat('SNR:',string(SNR(ExpsToPlot(iCond),28)));  
    end
    legend(LegendToBe)
    
    for iCond = 1:NrCond
        subplot(4,1,3)
        hold on
        plot(FTime{ExpsToPlot(iCond)}, FFiltDataAbs{ExpsToPlot(iCond)}(28,:)-VmTh{ExpsToPlot(iCond)},colorOrder{iCond})
        plot(FAHPxfit{ExpsToPlot(iCond),28}, FAHPFitted{ExpsToPlot(iCond),28}(:)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)}, colorOrder{iCond})
        plot(FAHPTime{ExpsToPlot(iCond)}(28), FAHPrmp{ExpsToPlot(iCond)}(28)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)},'rx')
        plot(FADPxfit{ExpsToPlot(iCond),28}, FADPFitted{ExpsToPlot(iCond),28}(:)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)}, colorOrder{iCond})
        plot(FADPTime{ExpsToPlot(iCond)}(28), FADPrmp{ExpsToPlot(iCond)}(28)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)},'ro')
        xlim([-3 6])
        ylim([-50 140])
        LegendToBe(iCond,:) = strcat('SNR:',string(SNR(ExpsToPlot(iCond),29)));  
    end
    legend(LegendToBe)
 
    for iCond = 1:NrCond
        subplot(4,1,4)
        hold on
        plot(FTime{ExpsToPlot(iCond)}, FFiltDataAbs{ExpsToPlot(iCond)}(29,:)-VmTh{ExpsToPlot(iCond)},colorOrder{iCond})
        plot(FAHPxfit{ExpsToPlot(iCond),29}, FAHPFitted{ExpsToPlot(iCond),29}(:)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)}, colorOrder{iCond})
        plot(FAHPTime{ExpsToPlot(iCond)}(29), FAHPrmp{ExpsToPlot(iCond)}(29)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)},'rx')
        plot(FADPxfit{ExpsToPlot(iCond),29}, FADPFitted{ExpsToPlot(iCond),29}(:)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)}, colorOrder{iCond})
        plot(FADPTime{ExpsToPlot(iCond)}(29), FADPrmp{ExpsToPlot(iCond)}(29)+VmRMP{ExpsToPlot(iCond)}-VmTh{ExpsToPlot(iCond)},'ro')
        xlim([-3 6])
        ylim([-50 140])
        LegendToBe(iCond,:) = strcat('SNR:',string(SNR(ExpsToPlot(iCond),30)));  
    end
    legend(LegendToBe)
 
    saveas(gcf,strcat(outputpath, 'Figures/MatlabExports/PairedAlignedTraces/',savename))
    pause
    close
    disp(strcat('Progress:  ' , string(iCell), '/', string(length(UniqueCells))', ' done'));
end

