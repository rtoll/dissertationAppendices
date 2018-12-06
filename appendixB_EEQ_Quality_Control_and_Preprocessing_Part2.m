clear all
close all
clc
startupEEGLAB
clear all
close all
clc
folderInput = 'F:\EEG_PREPROCESSING\STAGE1_RESAMPFILT';
folderOutput = 'F:\EEG_PREPROCESSING\STAGE2_BADCHANPAROXICA'
d = dir(folderInput);
filesInput = {d.name}';
clear d
filesInput = filesInput(endsWith(filesInput,'.mat'));
% pause for a small random amount of time to prevent overwriting when multiple instances of script by staggering
pause(rand()*10)
% load master cap definitions struct
load('D:\Dropbox\MATLAB\MatlabFunctions\matlab_master_definitions\masterEasyCapDefinitions.mat')
% define thresholds
threshFlatlineVoltage = .1; % uV
threshFlatlineTime = 1; % s
threshExtremaTime = 1; % s
threshExtremaMADmedian = 50; % uV
threshParoxRetainProp = .5; % proportion
threshParoxFreq = 5; % Hz
threshParoxMagTime = .5; % s
threshParoxStableTime = 1; % s
threshParoxSpikeT = 9; % T score
threshParoxStableT = 0; % T score
threshParoxInflation = .1; % T score
threshParoxMergeTime = 5; % s
threshCorrWindowTime = 1; % s
threshMagLimit = 100; % uV
for a=size(filesInput,1):-1:1
   fileInput = filesInput{a};
   fileOutput = fullfile(folderOutput,fileInput);
   % define a dummy placeholder to prevent overwriting when multiple instances of script running
   fileOutputDummy = strrep(fileOutput,'.mat','_DUMMY.mat');
   if ~exist(fileOutput,'file') && ~exist(fileOutputDummy,'file')
      disp('**************** PROCESSING ****************')
      disp(fileInput)
      disp('********************************************')
      % save dummy placeholder
      save(fileOutputDummy,'a')
      % load EEG
      load(fullfile(folderInput,fileInput))
      % identify and nullify flatline channels
      dataDiffMag = NaN(size(EEG.data));
      dataDiffMag(:,2:end) = abs(EEG.data(:,2:end) - EEG.data(:,1:end-1));
      dataDiffMagSubThresh = dataDiffMag < threshFlatlineVoltage;
      %
      chansFlatline = false(EEG.nbchan,1);
      for c=1:EEG.nbchan
         dataDiffMagSubThreshCC = bwconncomp(dataDiffMagSubThresh(c,:));
         if any(cellfun('length',dataDiffMagSubThreshCC.PixelIdxList) >= threshFlatlineTime * EEG.srate)
            chansFlatline(c) = true;
         end
      end
      EEG.data(chansFlatline,:) = NaN;
      % identify and nullify extrema channels
      dataMovMad = movmad(EEG.data,threshExtremaTime * EEG.srate,2);
      chansExtrema = nanmedian(dataMovMad,2) > threshExtremaMADmedian;
      EEG.data(chansExtrema,:) = NaN;
      % identify and remove paroxysmal data segments
      FILT = pop_eegfiltnew(EEG,threshParoxFreq,[]);
      data = FILT.data;
      clear FILT
      data(chansFlatline | chansExtrema,:) = NaN;
      dataMagSmooth = movmean(abs(data),threshParoxMagTime * EEG.srate,2);
      dataMagSmoothMedian = nanmedian(dataMagSmooth);
      pd = fitdist(dataMagSmoothMedian','gev');
      dataMagSmoothMedianGevNorm = (dataMagSmoothMedian - pd.mu) ./ pd.sigma;
      spikeIndices = find(dataMagSmoothMedianGevNorm > threshParoxSpikeT);
      paroxRetainProp = 0;
      k=0;
      while paroxRetainProp < threshParoxRetainProp
         parox = false(1,EEG.pnts);
         ccNotStable = bwconncomp(movmedian(dataMagSmoothMedianGevNorm,threshParoxStableTime * EEG.srate) > threshParoxStableT + (k * threshParoxInflation),4);
         for i=1:length(ccNotStable.PixelIdxList)
            segment = ccNotStable.PixelIdxList{i};
            if ~isempty(intersect(segment,spikeIndices))
               parox(segment) = true;
            end
         end
         parox = ~bwareaopen(~parox,threshParoxMergeTime * EEG.srate,4);
         paroxRetainProp = sum(~parox) / EEG.pnts;
         k=k+1;
         if k > 1000
            break
         end
      end
      EEG.parox = parox;
      EEG.paroxRetainProp = paroxRetainProp;
      ccRetain = bwconncomp(~parox,4);
      retainTimes = NaN(ccRetain.NumObjects,2);
      for i=1:ccRetain.NumObjects
         retainSegment = ccRetain.PixelIdxList{i};
         retainTimes(i,1) = retainSegment(1) / EEG.srate;
         retainTimes(i,2) = retainSegment(end) / EEG.srate;
      end
      RAW = EEG;
      RAW.data = EEG.rawdata;
      EEG = pop_select(EEG,'time',retainTimes);
      RAW = pop_select(RAW,'time',retainTimes);
      EEG.rawdata = RAW.data;
      clear RAW
      % identify channels with subthreshold maximal correlation
      dataWindowed = windowData(EEG.data,threshCorrWindowTime * EEG.srate);
      numWindows = size(dataWindowed,3);
      EEG.rhoWindow = NaN(EEG.nbchan,numWindows);
      for w=1:numWindows
         rho = abs(corrcoef(squeeze(dataWindowed(:,:,w))'));
         rho = diagNaN(rho);
         rho = sort(rho,2,'descend');
         rho = rho(:,~all(isnan(rho)));
         EEG.rhoWindow(:,w) = rho(:,1);
      end
      EEG.rhoWindowMedian = nanmedian(EEG.rhoWindow,2);
      if strcmp(EEG.paradigm,'REC')
         chansCorr = EEG.rhoWindowMedian < [master.thresh.rhoMedianREC];
      elseif strcmp(EEG.paradigm,'REO')
         chansCorr = EEG.rhoWindowMedian < [master.thresh.rhoMedianREO];
      end
      EEG.data(chansCorr,:) = NaN;
      % limit data maxima
      EEG.data = limitData(EEG.data,-threshMagLimit,threshMagLimit);
      % interpolate bad channels
      EEG.data = sphericalInterp(EEG.data,master.chanlocs);
      % rereference data
      EEG.meandata = mean(EEG.data);
      EEG.data = rereference(EEG.data,0,1);
      % assign channel categories
      chansBad = chansFlatline | chansExtrema | chansCorr;
      EEG.chans.bad = find(chansBad);
      EEG.chans.flatline = find(chansFlatline);
      EEG.chans.extrema = find(chansExtrema);
      EEG.chans.corr = find(chansCorr);
      EEG.chans.keep = find(~chansBad);
      % perform ICA, on GPU if available
      try
         gpu = gpuDevice(1);
         [weights,lrates,changes,angles,signs] = icaGPU(EEG.data(EEG.chans.keep,:),EEG.srate,gpu,true);
         EEG.ica.gpu = true;
      catch
         [weights,lrates,changes,angles,signs] = icaCPU(EEG.data(EEG.chans.keep,:),EEG.srate,true);
         EEG.ica.gpu = false;
      end
      EEG.ica.weights = weights;
      EEG.ica.lrates = lrates;
      EEG.ica.changes = changes;
      EEG.ica.angles = angles;
      EEG.ica.signs = signs;
      % assign preprocessing stage
      EEG.preprocessingStage = 2;
      % assign INFO variable for queries not needing to load data
      INFO = EEG;
      INFO = rmfield(INFO,'data');
      INFO = rmfield(INFO,'meandata');
      INFO = rmfield(INFO,'rawdata');
      INFO = rmfield(INFO,'times');
      % save EEG
      save(fileOutput,'EEG','INFO','-v7.3')
      % wait for system to register the save, when files are in the process of being written,
      % they have 0 bytes and are considered to not exist on some platforms
      while ~exist(fileOutput,'file')
      end
      % delete dummy placeholder
      delete(fileOutputDummy)
   end
   clearvars -except folder* files* master thresh* a
end