function [status] = porcupineStructuralRegistrationAndWarps(inputT1basename,basedir,fsldir)
   status = false;
   setenv('FSLDIR',fsldir);
   setenv('FSLOUTPUTTYPE', 'NIFTI');
   % define directories
   strawFolder = strcat(basedir,'/Processing/Stage2_StructuralRegistrationAndWarps/',inputT1basename,'_STRAW');
   % define filenames
   sffp = strcat(strawFolder,'/',inputT1basename,'_'); % strawFolder with filename prefix for strcat
   sdn = strcat(sffp,'SDN.nii'); % standardized and denoised
   natEB = strcat(sffp,'EB.nii'); % native extracted brain
   natEBM = strcat(sffp,'EBM.nii'); % native extracted brain mask
   mniEB = strcat(basedir,'/Processing/Templates/MNI/MNI_EXTRACTED_BRAIN_1MM.nii');
   mniEBM = strcat(basedir,'/Processing/Templates/MNI/MNI_EXTRACTED_BRAIN_MASK_1MM.nii');
   mni1mmT1 = strcat(basedir,'/Processing/Templates/MNI/MNI_T1_1MM.nii');
   mni2mmT1 = strcat(basedir,'/Processing/Templates/MNI/MNI_T1_2MM.nii');
   mniSCM = strcat(basedir,'/Processing/Templates/MNI/MNI_SUBCORT_MASK_1MM.nii');
   % initial registration
   system(sprintf('%s/flirt -in %s -out %s -ref %s -omat %s',strcat(fsldir,'/bin'),sdn,strcat(sffp,'firstFlirtTmpStage1'),mni1mmT1,strcat(sffp,'firstFlirtTmpStage1.mat')));
   % subcortical mask
   system(sprintf('%s/flirt -in %s -out %s -ref %s -omat %s -nosearch -refweight %s',strcat(fsldir,'/bin'),strcat(sffp,'firstFlirtTmpStage1'),strcat(sffp,'firstFlirtTmpStage2'),mni1mmT1,strcat(sffp,'firstFlirtTmpStage2.mat'),mniSCM));
   system(sprintf('%s/convert_xfm -omat %s -concat %s %s',strcat(fsldir,'/bin'),strcat(sffp,'firstFlirt.mat'),strcat(sffp,'firstFlirtTmpStage2.mat'),strcat(sffp,'firstFlirtTmpStage1.mat')));
   % generate registered images
   system(sprintf('%s/flirt -in %s -out %s -ref %s -applyxfm -init %s',strcat(fsldir,'/bin'),sdn,strcat(sffp,'firstFlirt'),mni1mmT1,strcat(sffp,'firstFlirt.mat')));
   % extracted brain mask
   system(sprintf('%s/flirt -in %s -out %s -ref %s -omat %s -nosearch -refweight %s',strcat(fsldir,'/bin'),strcat(sffp,'firstFlirtTmpStage1'),strcat(sffp,'firstFlirtTmpCortStage2'),mni1mmT1,strcat(sffp,'firstFlirtTmpCortStage2.mat'),mniEBM));
   system(sprintf('%s/convert_xfm -omat %s -concat %s %s',strcat(fsldir,'/bin'),strcat(sffp,'firstFlirtCort.mat'),strcat(sffp,'firstFlirtTmpCortStage2.mat'),strcat(sffp,'firstFlirtTmpStage1.mat')));
   system(sprintf('%s/flirt -in %s -out %s -ref %s -applyxfm -init %s',strcat(fsldir,'/bin'),sdn,strcat(sffp,'firstFlirtCort'),mni1mmT1,strcat(sffp,'firstFlirtCort.mat')));
   system(sprintf('%s/convert_xfm -omat %s -inverse %s',strcat(fsldir,'/bin'),strcat(sffp,'firstFlirtCortInv.mat'),strcat(sffp,'firstFlirtCort.mat')));
   system(sprintf('%s/flirt -in %s -ref %s -out %s -applyxfm -init %s',strcat(fsldir,'/bin'),mniEBM,sdn,natEBM,strcat(sffp,'firstFlirtCortInv.mat')));
   system(sprintf('%s/fslmaths %s -mas %s %s',strcat(fsldir,'/bin'),sdn,natEBM,natEB));
   % execute flirt registration sequence
   system(sprintf('%s/flirt -ref %s -in %s -out %s -omat %s',strcat(fsldir,'/bin'),mniEB,natEB,strcat(sffp,'nativeToMni'),strcat(sffp,'nativeToMni.mat')));
   system(sprintf('%s/convert_xfm -omat %s -inverse %s',strcat(fsldir,'/bin'),strcat(sffp,'mniToNative.mat'),strcat(sffp,'nativeToMni.mat')));
   % execute level 1 fnirt registration sequence
   system(sprintf('%s/fnirt --ref=%s --in=%s --aff=%s --cout=%s --iout=%s --intout=%s --subsamp=4,4,2,2,1,1 --miter=5,5,5,5,5,10 --infwhm=8,6,5,4.5,3,2 --reffwhm=8,6,5,4,2,0 --lambda=300,150,100,50,40,30 --estint=1,1,1,1,1,0 --applyrefmask=1,1,1,1,1,1 --applyinmask=1 --intmod=global_non_linear_with_bias',strcat(fsldir,'/bin'),mni2mmT1,sdn,strcat(sffp,'nativeToMni.mat'),strcat(sffp,'nativeToMniWarp1'),strcat(sffp,'nativeToMniWarped1'),strcat(sffp,'nativeToMniIntensities')));
   % execute level 2 fnirt registration sequence
   system(sprintf('%s/fnirt --ref=%s --in=%s --inwarp=%s --intin=%s --cout=%s --iout=%s --subsamp=1 --miter=20 --infwhm=2 --reffwhm=0 --lambda=30 --estint=0 --warpres=6,6,6 --intmod=global_non_linear_with_bias --applyrefmask=1 --applyinmask=1 --minmet=scg --jacrange=-5,100',strcat(fsldir,'/bin'),mni2mmT1,sdn,strcat(sffp,'nativeToMniWarp1'),strcat(sffp,'nativeToMniIntensities'),strcat(sffp,'nativeToMniWarp2'),strcat(sffp,'nativeToMniWarped2')));
   % execute level 3 fnirt registration sequence
   system(sprintf('%s/fnirt --ref=%s --in=%s --inwarp=%s --intin=%s --cout=%s --iout=%s --subsamp=1 --miter=15 --infwhm=2 --reffwhm=0 --lambda=10 --estint=0 --warpres=2,2,2 --intmod=global_non_linear_with_bias --applyrefmask=1 --applyinmask=1 --minmet=scg --jacrange=-5,100',strcat(fsldir,'/bin'),mni2mmT1,sdn,strcat(sffp,'nativeToMniWarp2'),strcat(sffp,'nativeToMniIntensities'),strcat(sffp,'nativeToMniWarp'),strcat(sffp,'nativeToMniWarped')));
   % apply warp fields to native extracted brain
   system(sprintf('%s/invwarp -w %s -o %s -r %s',strcat(fsldir,'/bin'),strcat(sffp,'nativeToMniWarp'),strcat(sffp,'mniToNativeWarp'),natEB));
   system(sprintf('%s/convert_xfm -omat %s -inverse %s',strcat(fsldir,'/bin'),strcat(sffp,'mniToNative.mat'),strcat(sffp,'nativeToMni.mat')));
   system(sprintf('%s/applywarp -w %s -i %s -r %s -o %s -d float',strcat(fsldir,'/bin'),strcat(sffp,'mniToNativeWarp'),mniEBM,sdn,strcat(sffp,'brainFnirtMask')));
   system(sprintf('%s/fslmaths %s -thr 0 -bin %s -odt short',strcat(fsldir,'/bin'),strcat(sffp,'brainFnirtMask'),strcat(sffp,'brainFnirtMask')));
   system(sprintf('%s/fslmaths %s -mas %s %s',strcat(fsldir,'/bin'),sdn,strcat(sffp,'brainFnirtMask'),strcat(sffp,'brainFnirt')));
   % delete unnecessary files
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtTmpStage1.mat')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtTmpStage1.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtTmpStage2.mat')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtTmpStage2.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirt.mat')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirt.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtTmpCortStage2.mat')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtCort.mat')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtTmpCortStage2.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtCort.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'firstFlirtCortInv.mat')));
   system(sprintf('rm -f %s',strcat(sffp,'nativeToMniWarp1.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'nativeToMniWarped1.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'nativeToMniWarp2.nii')));
   system(sprintf('rm -f %s',strcat(sffp,'nativeToMniWarped2.nii')));
   status = true;
end