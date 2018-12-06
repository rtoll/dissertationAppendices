function [isValidFilename] = porcupineValidateFilenameConformsWithLabStandard(inputT1)
   isValidFilename = true;
   try
      lastSlash = strfind(inputT1,'/');
      if isempty(lastSlash)
         lastSlash = 0;
      else
         lastSlash = lastSlash(end);
      end
      niiExt = strfind(inputT1,'.nii');
      % make sure extension is only .nii and not .nii.gz
      if ~isempty(strfind(inputT1,'.nii.gz'))
         isValidFilename = false;
      end
      baseName = inputT1(lastSlash + 1:niiExt - 1);
      delims = strfind(baseName,'_');
      % there must be exactly 3 underscores
      if length(delims) ~= 3
         isValidFilename = false;
      end
      studyName = baseName(1:delims(1) - 1);
      % studyName must be uppercase
      if ~strcmp(studyName,upper(studyName))
         isValidFilename = false;
      end
      studyID = baseName(delims(1) + 1:delims(2) - 1);
      % studyID must be 3 or more numbers long, use leading zeros if necessary
      if length(studyID) < 3
         isValidFilename = false;
      end
      % studyID numbers must be numbers
      if isnan(str2double(studyID))
         isValidFilename = false;
      end
      redcapID = baseName(delims(2) + 1:delims(3) - 1);
      % recapID must be preceded by RC
      if ~strcmp(redcapID(1:2),'RC')
         isValidFilename = false;
      end
      redcapID = redcapID(3:end);
      % redcapID numbers must be exactly 4 numbers long
      if length(redcapID) ~= 4
         isValidFilename = false;
      end
      % redcapID numbers must be numbers
      if isnan(str2double(redcapID))
         isValidFilename = false;
      end
      % baseName must end with T1
      if ~strcmp(baseName(end - 1:end),'T1')
         isValidFilename = false;
      end
   catch
      isValidFilename = false;
   end
end