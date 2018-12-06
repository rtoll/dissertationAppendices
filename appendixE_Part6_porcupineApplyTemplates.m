function [status] = porcupineApplyTemplates(inputT1basename,basedir,fsldir)
   status = false;
   setenv('FSLDIR',fsldir);
   setenv('FSLOUTPUTTYPE', 'NIFTI');
   % define and create directories
   strawFolder = strcat(basedir,'/Processing/Stage2_StructuralRegistrationAndWarps/',inputT1basename,'_STRAW');
   ntFolder = strcat(basedir,'/Processing/Stage3_ApplyTemplates/',inputT1basename,'_TEMPLATES'); % native templates folder
   ttFolder = strcat(basedir,'/Processing/Templates/Targeting'); % targeting templates folder
   system(sprintf('mkdir -p %s',ntFolder));
   % define filenames
   sffp = strcat(strawFolder,'/',inputT1basename,'_'); % strawFolder with filename prefix for strcat
   sdn = strcat(sffp,'SDN'); % standardized and denoised
   mniToNativeWarp = strcat(sffp,'mniToNativeWarp');
   ntffp = strcat(ntFolder,'/',inputT1basename,'_'); % ntFolder with filename prefix for strcat
   sdnnzb = strcat(ntffp,'SDNNZB'); % standardized and denoised T1 with non-zero background
   % add 1 to standardized and denoised T1 to prevent null warping
   system(sprintf('%s/fslmaths %s -add 1 %s -odt short',strcat(fsldir,'/bin'),sdn,sdnnzb));
   % obtain the x dimension of native reference for rad to neuro coord conversion
   [~,dim1] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',sdn,'dim1'));
   dim1 = str2double(dim1);
   % assemble the list of templates to process
   templates = dir(ttFolder);
   templates = {templates.name};
   keepTemplates = {};
   j=1;
   for i=1:size(templates,2);
      template = templates{i};
      if ~isempty(strfind(template,'.nii'))
         keepTemplates(j) = {template(1:strfind(template,'.nii') - 1)};
         j=j+1;
      end
   end
   % apply warp fields to mni templates and write target matrix
   system(sprintf('export FSLOUTPUTTYPE=NIFTI'));
   templates = keepTemplates;
   j=1;
   for templateIndex = 1:size(templates,2)
      template = strcat(ttFolder,'/',templates{templateIndex});
      clc
      fprintf('warping MNI templates %u%% complete\n',round(100 * templateIndex / size(templates,2)))
      templateInNativeSpace = strcat(ntffp,templates{templateIndex});
      system(sprintf('%s/applywarp -i %s -r %s -w %s -o %s',strcat(fsldir,'/bin'),template,sdnnzb,mniToNativeWarp,templateInNativeSpace));
      targNum = str2double(template(strfind(template,'_') - 3 : strfind(template,'_') - 1));
      [~,maxVox] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslstats',templateInNativeSpace,'-x'));
      % delete template when complete
      system(sprintf('rm -f %s.nii',templateInNativeSpace));
      maxVox = str2num(maxVox);
      vtmRad(j,:) = [targNum,maxVox]; % vtmRad is in radiological RLPAISvoxel format
      maxVox = [dim1 - maxVox(1) + 1,maxVox(2),maxVox(3)];
      vtm(j,:) = [targNum,maxVox]; % vtm is in neuro LRPAISvoxel format
      j=j+1;
   end
   save(strcat(ntffp,'VTMrad'),'vtmRad') % voxel targets matrix
   save(strcat(ntffp,'VTM'),'vtm') % voxel targets matrix
   status = true;
end