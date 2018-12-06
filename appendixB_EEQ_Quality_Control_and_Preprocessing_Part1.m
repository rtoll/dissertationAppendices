clear all
close all
clc
folderInput = 'F:\EEG_PREPROCESSING\STAGE0_RAW';
folderOutput = 'F:\EEG_PREPROCESSING\STAGE1_RESAMPFILT';
d = dir(folderInput);
filesInput = {d.name}';
clear d
filesInput = filesInput(endsWith(filesInput,'.mat'));
% pause for a small random amount of time to prevent overwriting when multiple instances of script by staggering
pause(rand()*10)
for a=1:size(filesInput,1)
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
      % determine and assign studyName
      EEG.studyName = fileInput(1 : underscoreIndex(fileInput,1) - 1);
      % determine and assign subjectID
      EEG.subjectID = fileInput(underscoreIndex(fileInput,1) + 1 : underscoreIndex(fileInput,2) - 1);
      % determine and assign redcapID
      rcIndex = strfind(fileInput,'_RC');
      if ~isempty(rcIndex)
         EEG.redcapID = fileInput(rcIndex + 1 : rcIndex + 1 + 5);
      else
         EEG.redcapID = '';
      end
      % determine and assign visitNum
      visitIndex = strfind(fileInput,'_VISIT');
      if ~isempty(visitIndex)
         EEG.visit = fileInput(visitIndex + 6);
      else
         EEG.visit = '';
      end
      % determine and assign partNum
      partIndex = strfind(fileInput,'_PART');
      if ~isempty(partIndex)
         EEG.part = fileInput(partIndex + 5);
      else
         EEG.part = '';
      end
      % determine and assign paradigm
      if contains(fileInput,'_REC')
         EEG.paradigm = 'REC';
      elseif contains(fileInput,'_REO')
         EEG.paradigm = 'REO';
      else
         error('bad paradigm')
      end
      % determine and assign dateVec
      EEG.dateVec = datevec(EEG.datetimestamp,'yyyymmddHHMM');
      % determine and assign dateNum
      EEG.dateNum = datenum(EEG.dateVec);
      % determine and assign impedances
      EEG.imps = NaN(EEG.nbchan + 2,1);
      foundImpedance = find(startsWith(EEG.vhdr,'Impedance'),1,'last');
      if ~isempty(foundImpedance)
         vhdrImps = EEG.vhdr(foundImpedance + 1 : size(EEG.vhdr,1));
         vhdrImpsLabels = extractBefore(vhdrImps,':');
         vhdrImpsValues = str2double(extractAfter(vhdrImps,':'));
         vhdrImpsValues(isempty(vhdrImpsValues)) = NaN;
         if length(vhdrImpsValues) > EEG.nbchan
            EEG.imps(1 : EEG.nbchan) = vhdrImpsValues(1 : EEG.nbchan);
         end
         if any(strcmp(vhdrImpsLabels,'Ref'))
            EEG.imps(EEG.nbchan + 1) = vhdrImpsValues(strcmp(vhdrImpsLabels,'Ref'));
         end
         if any(strcmp(vhdrImpsLabels,'Gnd'))
            EEG.imps(EEG.nbchan + 2) = vhdrImpsValues(strcmp(vhdrImpsLabels,'Gnd'));
         end
      end
      % ensure data is class double precision
      if ~isa(EEG.data,'double')
         EEG.data = double(EEG.data);
      end
      % clear events
      EEG.event = [];
      % resample
      if EEG.srate > 256
         EEG = pop_resample(EEG,250);
      end
      % ensure data is class double precision again, as certain EEGLAB functions revert to single
      if ~isa(EEG.data,'double')
         EEG.data = double(EEG.data);
      end
      % save raw data in temporary RAW struct
      RAW = EEG;
      % notch filter line noise
      EEG = pop_eegfiltnew(EEG,58,62,[],1);
      if ~isa(EEG.data,'double')
         EEG.data = double(EEG.data);
      end
      % notch filter line noise harmonic
      EEG = pop_eegfiltnew(EEG,118,122,[],1);
      % ensure data is class double precision again, as certain EEGLAB functions revert to single
      if ~isa(EEG.data,'double')
         EEG.data = double(EEG.data);
      end
      % high pass filter
      EEG = pop_eegfiltnew(EEG,2,[]);
      % ensure data is class double precision again, as certain EEGLAB functions revert to single
      if ~isa(EEG.data,'double')
         EEG.data = double(EEG.data);
      end
      % low pass filter
      EEG = pop_eegfiltnew(EEG,[],48);
      % ensure data is class double precision again, as certain EEGLAB functions revert to single
      if ~isa(EEG.data,'double')
         EEG.data = double(EEG.data);
      end
      % trim edges of data to remove filtering artifacts
      EEG = pop_select(EEG,'time',[EEG.xmin + 2 , EEG.xmax - 2]);
      RAW = pop_select(RAW,'time',[RAW.xmin + 2 , RAW.xmax - 2]);
      % ensure data is class double precision again, as certain EEGLAB functions revert to single
      if ~isa(EEG.data,'double') || ~isa(RAW.data,'double')
         EEG.data = double(EEG.data);
         RAW.data = double(RAW.data);
      end
      % save raw data in EEG struct field
      EEG.rawdata = RAW.data;
      clear RAW
      % assign preprocessing stage
      EEG.preprocessingStage = 1;
      % assign INFO variable for queries not needing to load data
      INFO = EEG;
      INFO = rmfield(INFO,'data');
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
   clearvars -except folder* files* a
end