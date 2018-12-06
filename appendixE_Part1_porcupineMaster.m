function [status] = porcupineMaster(fsldir,basedir)
% THIS IS THE MASTER CALL FUNCTION FOR THE PORCUPINE PIPELINE
   % add the admin path
   status = false;
   addpath(genpath(basedir))
   system(sprintf('chmod -R 777 %s',strcat(basedir,'/InputT1/')));
   % preamble
   porcupineArt(basedir);
   pause(1)
   % assemble the list of T1s to process
   inputT1s = dir(strcat(basedir,'/InputT1'));
   inputT1s = {inputT1s.name};
   keepInputT1s = {};
   j=1;
   for i=1:size(inputT1s,2)
      inputT1 = inputT1s{i};
      if ~strcmp(inputT1(1),'.')
         keepInputT1s(j) = strcat(basedir,'/InputT1/',inputT1s(i));
         j=j+1;
      end
   end
   if ~isempty(keepInputT1s)
      inputT1 = keepInputT1s{1};
      repeat = true;
   else
      repeat = false;
      disp('DID NOT FIND ANY T1s TO PROCESS')
      pause(5)
   end
   while repeat
      % proceed through targeting pipeline
      repeat = false;
      %  ensure filename is valid
      disp('PHASE 1 of 8: VALIDATE FILENAME')
      isValidFilename = porcupineValidateFilenameConformsWithLabStandard(inputT1);
      if isValidFilename         
         % reduce inputT1 to its basename
         inputT1basename = inputT1(length(strcat(basedir,'/InputT1/')) + 1 : strfind(inputT1,'RC') + 5);
         moveFailedSubjectToQuarantine = false;
         clc
         disp(strcat('PROCESSING: ',inputT1basename))
         disp('PHASE 1 of 8 COMPLETE')
         pause(1)
         clc            
         % standardize the format of inputT1 to comply with visor requirements
         if ~moveFailedSubjectToQuarantine
            disp('PHASE 2 of 8: STANDARDIZE T1 TO VISOR COMPATIBLE FORMAT')
            disp(inputT1basename)
            try
               porcupineStandardizeT1(inputT1basename,basedir,fsldir);
               disp('PHASE 2 of 8 COMPLETE')
               pause(1)
               clc
            catch
               disp('PHASE 2 of 8 FAILED, MOVING SUBJECT TO QUARANTINE')
               moveFailedSubjectToQuarantine = true;
               failPhase = 2;
               pause(1)
               clc
            end
         end            
         % denoise resulting T1 which has suffix '_OVF' (original in visor format)
         if ~moveFailedSubjectToQuarantine
            disp('PHASE 3 of 8: DENOISE T1')
            disp(inputT1basename)
            try
               porcupineDenoiseT1(inputT1basename,basedir,fsldir);
               disp('PHASE 3 of 8 COMPLETE')
               pause(1)
               clc
            catch
               disp('PHASE 3 of 8 FAILED, MOVING SUBJECT TO QUARANTINE')
               moveFailedSubjectToQuarantine = true;
               failPhase = 3;
               pause(1)
               clc
            end
         end
         % execute structural processing sequence            
         if ~moveFailedSubjectToQuarantine
            disp('PHASE 4 of 8: STRUCTURAL REGISTRATION AND WARP COMPUTATION')
            disp(inputT1basename)
            try
               porcupineStructuralRegistrationAndWarps(inputT1basename,basedir,fsldir);
               disp('PHASE 4 of 8 COMPLETE')
               pause(1)
               clc
            catch
               disp('PHASE 4 of 8 FAILED, MOVING SUBJECT TO QUARANTINE')
               moveFailedSubjectToQuarantine = true;
               failPhase = 4;
               pause(1)
               clc
            end
         end
         % apply warp fields to mni templates and write target matrix            
         if ~moveFailedSubjectToQuarantine
            disp('PHASE 5 of 8: WARP MNI TARGET TEMPLATES INTO NATIVE SPACE')
            disp(inputT1basename)
            try
               porcupineApplyTemplates(inputT1basename,basedir,fsldir);
               disp('PHASE 5 of 8 COMPLETE')
               pause(1)
               clc
            catch
               disp('PHASE 5 of 8 FAILED, MOVING SUBJECT TO QUARANTINE')
               moveFailedSubjectToQuarantine = true;
               failPhase = 5;
               pause(1)
               clc
            end
         end
         % refine initial fiducials estimate to actual positions            
         if ~moveFailedSubjectToQuarantine
            disp('PHASE 6 of 8: REFINE FIDUCIAL MARKER LOCATIONS')
            disp(inputT1basename)
            try
               porcupineRefineFiducials(inputT1basename,basedir,fsldir);
               disp('PHASE 6 of 8 COMPLETE')
               pause(1)
               clc
            catch
               disp('PHASE 6 of 8 FAILED, MOVING SUBJECT TO QUARANTINE')
               moveFailedSubjectToQuarantine = true;
               failPhase = 6;
               pause(1)
               clc
            end
         end
         % calculate electrode locations            
         if ~moveFailedSubjectToQuarantine
            disp('PHASE 7 of 8: CALCULATE ELECTRODE POSITIONS AND WRITE IMAGE AND CHANNEL DATA')
            disp(inputT1basename)
            try
               porcupineApplyElectrodes(inputT1basename,basedir,fsldir);
               disp('PHASE 7 of 8 COMPLETE')
               pause(1)
               clc
            catch
               disp('PHASE 7 of 8 FAILED, MOVING SUBJECT TO QUARANTINE')
               moveFailedSubjectToQuarantine = true;
               failPhase = 7;
               pause(1)
               clc
            end
         end
         % clean up, remove temporary files, assemble outputs and move to output folder
         if ~moveFailedSubjectToQuarantine
            disp('PHASE 8 of 8: REMOVE TEMPORARY FILES AND ASSEMBLE FINAL OUTPUT')
            porcupineFinalCleanUp(inputT1basename,basedir);   
            disp('PHASE 8 of 8 COMPLETE')
            disp(sprintf('%s PORCUPINE PIPELINE COMPLETE',inputT1basename))
            pause(2)
         else
            % define and create directory
            quarantineFolder = strcat(basedir,'/Quarantine/',inputT1basename);
            system(sprintf('mkdir -p %s',quarantineFolder));
            % define filename prefix
            qfp = strcat(quarantineFolder,'/',inputT1basename,'_');
            fclose('all');
            fid = fopen(strcat(qfp,'FAILED_AT_PHASE_',num2str(failPhase),'.txt'),'w');
            fwrite(fid,strcat(qfp,'FAILED_AT_PHASE_',num2str(failPhase)));
            fclose('all');
            % move inputT1
            system(sprintf('mv -f %s %s.nii',inputT1,strcat(qfp,'RAW')));
            % move files and folders to quarantine if they exist
            if exist(strcat(basedir,'/Processing/Stage1_StandardizeT1/',inputT1basename,'_OVF.nii'),'file')
               system(sprintf('mv -f %s %s.nii',strcat(basedir,'/Processing/Stage1_StandardizeT1/',inputT1basename,'_OVF.nii'),strcat(qfp,'OVF')));
            end
            if exist(strcat(basedir,'/Processing/Stage2_StructuralRegistrationAndWarps/',inputT1basename,'_STRAW'),'dir')
               system(sprintf('mv -f %s %s',strcat(basedir,'/Processing/Stage2_StructuralRegistrationAndWarps/',inputT1basename,'_STRAW'),strcat(qfp,'STRAW')))
            end
            if exist(strcat(basedir,'/Processing/Stage3_ApplyTemplates/',inputT1basename,'_TEMPLATES'),'dir')
               system(sprintf('mv -f %s %s',strcat(basedir,'/Processing/Stage3_ApplyTemplates/',inputT1basename,'_TEMPLATES'),strcat(qfp,'TEMPLATES')))
            end
         end
      else
         % if invalid filename, display warning and continue through the remainder of the queue
         disp(sprintf('WARNING: %s does NOT comply with lab naming convention, refer to HOW_TO_USE_PORCUPINE_READ_ME.txt\ncontinuing with the next T1 in queue',inputT1))
         pause(10)
         clc
      end
      % assemble the list of remaining T1s to process
      inputT1s = dir(strcat(basedir,'/InputT1'));
      inputT1s = {inputT1s.name};
      keepInputT1s = {};
      j=1;
      for i=1:size(inputT1s,2);
         inputT1 = inputT1s{i};
         if ~strcmp(inputT1(1),'.')
            keepInputT1s(j) = strcat(basedir,'/InputT1/',inputT1s(i));
            j=j+1;
         end
      end
      if ~isempty(keepInputT1s)
         inputT1 = keepInputT1s{1};
         repeat = true;
      else
         disp('DID NOT FIND ADDITIONAL T1s TO PROCESS, THE PORCUPINE IS GOING BACK TO ITS DEN NOW')
         pause(2)
      end
   end
   status = true;
   exit
end