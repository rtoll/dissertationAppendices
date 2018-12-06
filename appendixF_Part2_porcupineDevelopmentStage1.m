clear all
close all
clc

baseDir = '/Users/russelltoll/Dropbox/hedgehog/PORCUPINE';
fslDirectory = '/usr/local/fsl/bin';

% remove above in actual call stack

%% Stage 1. Validate filenames of T1s in InputT1 directory
stage = '1.0';

% Stage 1.1. 
stage = '1.1';
inputT1s = dir(fullfile(baseDir,'InputT1'));
inputT1s = {inputT1s.name}';
inputT1s = inputT1s(contains(inputT1s,'.nii') & ~contains(inputT1s,'INVALID'));
isValidFilename = false(size(inputT1s));
% Stage 1.2.
stage = '1.2';
for i=1:size(inputT1s,1)
   inputT1 = inputT1s{i};
   try
      % Stage 1.2a
      assert(length(strfind(inputT1,'_')) == 3 || length(strfind(inputT1,'_')) == 4,'stage1:step2a','there must be either 3 (standard) or 4 (multi visit) underscores in the filename')
      % Stage 1.2b
      studyName = inputT1(1:underscoreIndex(inputT1,1)-1);
      assert(isequal(1:length(studyName),regexp(studyName,'[A-Z]')),'stage1:step2b','study name must contain only upper case letters')
      % Stage 1.2c
      subjectNumber = inputT1(underscoreIndex(inputT1,1)+1:underscoreIndex(inputT1,2)-1);
      assert(length(subjectNumber) >= 3 && isequal(1:length(subjectNumber),regexp(subjectNumber,'[0-9]')),'stage1:step2c','subject number must contain only digits and must be at least 3 digits long (use leading zeros if necessary)')
      % Stage 1.2d
      redcapNumber = inputT1(underscoreIndex(inputT1,2)+1:underscoreIndex(inputT1,3)-1);
      assert(regexp(redcapNumber,'RC[0-9][0-9][0-9][0-9]') == 1,'stage1:step2d','redcap number must have format RC####')
      % Stage 1.2e
      if length(strfind(inputT1,'_')) == 4
         visitNumber = inputT1(underscoreIndex(inputT1,3)+1:underscoreIndex(inputT1,4)-1);
         assert(regexp(visitNumber,'VISIT[1-9]') == 1,'stage1:step2e','visit identifier must have format VISIT#')
      end
      % Stage 1.2f
      assert(endsWith(inputT1,'_T1.nii'),'stage1:step2f','filename must end with _T1.nii (not nii.gz)')
      isValidFilename(i) = true;
      
      % add check only positive
      
   catch errorReportStage1
      disp('******* WARNING *******')
      disp('******* INVALID FILENAME, SKIPPING *******')
      disp(inputT1)
      disp(errorReportStage1.identifier)
      disp(errorReportStage1.message)
   end
end
invalidT1s = inputT1s(~isValidFilename);
for i=1:size(invalidT1s,1)   
   movefile(fullfile(baseDir,'InputT1',invalidT1s{i}),fullfile(baseDir,'InputT1',strcat('INVALID_NAME_',invalidT1s{i})));
end
inputT1s = inputT1s(isValidFilename);

% Stage 1 cleanup
clearvars -except baseDir logFilename fslDirectory inputT1s
