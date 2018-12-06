function [status] = porcupineStandardizeT1(inputT1basename,basedir,fsldir)
   % assumes inputT1 has passed valid filename check
   status = false;
   setenv('FSLDIR',fsldir);
   setenv('FSLOUTPUTTYPE', 'NIFTI');
   % define directories
   rawFolder = strcat(basedir,'/Processing/Stage0_RawT1');
   stanFolder = strcat(basedir,'/Processing/Stage1_StandardizeT1');
   % define filenames
   rfp = strcat(rawFolder,'/',inputT1basename,'_'); % rawFolder with filename prefix for strcat
   sfp = strcat(stanFolder,'/',inputT1basename,'_'); % stanFolder with filename prefix for strcat
   tempNoNans = strcat(sfp,'tempNoNans');
   tempReorient = strcat(sfp,'tempReorient');
   tempRadiological = strcat(sfp,'tempRadiological');
   tempInt16 = strcat(sfp,'tempInt16');
   ovf = strcat(sfp,'OVF'); % OVF = original T1 in visor format
   % move inputT1 to rawT1 folder so the processing queue can be updated for multiple nodes
   system(sprintf('mv -f %s.nii %s.nii',strcat(basedir,'/InputT1/',inputT1basename,'_T1'),strcat(rfp,'T1')));
   % convert NaNs (if present) to zero
   system(sprintf('%s/%s %s %s %s',strcat(fsldir,'/bin'),'fslmaths',strcat(rfp,'T1'),'-nan',tempNoNans));
   % reorient to standard orientation (RL PA IS)
   system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslreorient2std',tempNoNans,tempReorient));
   % ensure radiological format
   [~,hdrOrient] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslorient','-getorient',tempReorient));
   if strcmp(strtrim(hdrOrient),'NEUROLOGICAL')
      system(sprintf('%s/%s %s %s %s',strcat(fsldir,'/bin'),'fslswapdim',tempReorient,'-x y z',tempRadiological));
      system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslorient','-forceradiological',tempRadiological));
   else
      system(sprintf('cp -f %s.nii %s.nii',tempReorient,tempRadiological));
   end
   % convert to data type int16
   system(sprintf('%s/%s %s %s -odt short',strcat(fsldir,'/bin'),'fslmaths',tempRadiological,tempInt16));
   % make ovf in preparation for denoising
   system(sprintf('cp -f %s.nii %s.nii',tempInt16,ovf));
   % remove temporary niftis
   system(sprintf('rm -f %s.nii',tempNoNans));
   system(sprintf('rm -f %s.nii',tempReorient));
   system(sprintf('rm -f %s.nii',tempRadiological));
   system(sprintf('rm -f %s.nii',tempInt16));
   status = true;
end