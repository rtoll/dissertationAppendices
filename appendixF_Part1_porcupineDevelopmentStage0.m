clear all
close all
clc

baseDir = '/scratch/users/rtoll/PORCUPINE/';

% remove above in actual call stack

% try % main

% Stage 0. Preamble and procedural checks
stage = '0.0';

% Stage 0.1. Initialize log
stage = '0.1';
timestamp = datevec(now);
logFilename = fullfile(baseDir,'Logs',sprintf('porcupineLog_%04u_%02u_%02u_%02u_%02u_%02u_%03u.log',timestamp(1),timestamp(2),timestamp(3),timestamp(4),timestamp(5),floor(timestamp(6)),round(1000*rem(timestamp(6),floor(timestamp(6))))));
diary(logFilename)

% Stage 0.2. Display title screen
stage = '0.2';
% display porc
disp('Stage 0. Preamble and procedural checks')

% Stage 0.3. Check operating system is unix
stage = '0.3';
assert(isunix,'stage0:step3','operating system must be unix')

% Stage 0.4. Log user and machine information
stage = '0.4';
[~,userInfo] = system('id');
disp(userInfo)
[~,machineInfo] = system('hostname');
disp(machineInfo)

% Stage 0.5. Check Matlab version and toolboxes
stage = '0.5';
matlabVersion = ver;
matlabYear = matlabVersion(strcmp({matlabVersion.Name},'MATLAB')).Release;
matlabYear = str2double(matlabYear(3:6));
assert(matlabYear >= 2017,'stage0:step5a','this program requires MATLAB (2017a) or later')
assert(any(strcmp({matlabVersion.Name},'Image Processing Toolbox')),'stage0:step5b','this program requires the Image Processing Toolbox')

% Stage 0.6. Check FSL version and availability
stage = '0.6';
[terminalCheck,fslDirectory] = system('echo $FSLDIR');
assert(terminalCheck == 0,'stage0:step6a','problem with terminal returning cmdout')
try % fsl accesible
   assert(contains(fslDirectory,'fsl'),'stage0:step6a','FSL is not accessible to Matlab, check module load protocols and shell')
catch
   if exist('/usr/local/fsl','dir')
      fslDirectory = '/usr/local/fsl';
   else
      assert(contains(fslDirectory,'fsl'),'stage0:step6a','FSL is not accessible to Matlab, check module load protocols and shell')
   end
end % fsl accesible
fslVersion = fslDirectory(slashIndex(fslDirectory,-1)+1:end)
assert(fslVersion >= 5,'stage0:step6b','this program requires FSL version 5 or later')
% append /bin to fslDirectory for calls to FSL functions
fslDirectory = strcat(fslDirectory,'/bin');

% Stage 0.7. Check all requisite directories exist, read me file is present, and unauthorized files/folders are not present (PORCUPINE folder has not been altered)
stage = '0.7';
dirInfo = dir(baseDir);
authorizedFolders = {'.','..','Admin','InputT1','Logs','Output','Processing','Quarantine'};
foldersPresent = {dirInfo([dirInfo.isdir]).name};
filesPresent = {dirInfo(~[dirInfo.isdir]).name};
unauthorizedFolders = setdiff(foldersPresent,authorizedFolders);
unauthorizedFiles = setdiff(filesPresent,'PORCUPINE_READ_ME.txt');
assert(isempty(setdiff(authorizedFolders,foldersPresent)),'stage0:step7a','a requisite subfolder is missing from the PORCUPINE directory')
assert(exist(fullfile(baseDir,'PORCUPINE_READ_ME.txt'),'file')==2,'stage0:step7b','the PORCUPINE_READ_ME file is required to be present')
if sum(~[dirInfo.isdir]) ~= 1 || isempty(setdiff(foldersPresent,authorizedFolders))
   warning('additional files or folders are not authorized in the PORCUPINE directory, refer to the PORCUPINE_READ_ME.txt file. Removing them.')
   mkdir(strcat(baseDir,'/Quarantine/UnauthorizedFilesAndFolders'))
   for i=1:size(unauthorizedFiles,2)
      unauthorizedFile = fullfile(baseDir,unauthorizedFiles{i});
      movefile(unauthorizedFile,strrep(unauthorizedFile,baseDir,strcat(baseDir,'/Quarantine/UnauthorizedFilesAndFolders')));
   end
   for i=1:size(unauthorizedFolders,2)
      unauthorizedFolder = fullfile(baseDir,unauthorizedFolders{i});
      movefile(unauthorizedFolder,strrep(unauthorizedFolder,baseDir,strcat(baseDir,'/Quarantine/UnauthorizedFilesAndFolders')));
   end
   rmdir(strcat(baseDir,'/Quarantine/UnauthorizedFilesAndFolders'),'s')
end
% if removal of unauthorized files and/or folders was unsuccessful, throw error
assert(sum(~[dirInfo.isdir]) == 1,'stage0:step7c','no other files are authorized in the PORCUPINE directory')
assert(isempty(setdiff(foldersPresent,authorizedFolders)),'stage0:step7d','no additional subfolders are authorized in the PORCUPINE directory')

% Stage 0 cleanup
clearvars -except baseDir logFilename fslDirectory

diary off
movefile(logFilename,strrep(logFilename,'.log','_COMPLETE.log'));
disp('program complete')
disp('terminating program')
% exit

%%
% catch errorReportMain % main
% 
% disp('******* ERROR *******')
% fprintf('error in stage %s\n',stage)
% disp(errorReportMain.identifier)
% disp(errorReportMain.message)
% diary off
% movefile(logFilename,strrep(logFilename,'.log','_ERROR.log'));
% disp('terminating program')
% % exit
%    
% end % main






