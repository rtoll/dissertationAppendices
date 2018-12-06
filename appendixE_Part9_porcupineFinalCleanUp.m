function [status] = porcupineFinalCleanUp(inputT1basename,basedir)
   status = false;
   % define and create directories
   rawFolder = strcat(basedir,'/Processing/Stage0_RawT1');
   strawFolder = strcat(basedir,'/Processing/Stage2_StructuralRegistrationAndWarps/',inputT1basename,'_STRAW');
   ntFolder = strcat(basedir,'/Processing/Stage3_ApplyTemplates/',inputT1basename,'_TEMPLATES');
   outFolder = strcat(basedir,'/Output/',inputT1basename,'_PORCUPINE');
   system(sprintf('mkdir -p %s',outFolder));
   outSTRAWfolder = strcat(basedir,'/Output/',inputT1basename,'_PORCUPINE/',inputT1basename,'_STRAW');
   outICDfolder = strcat(basedir,'/Output/',inputT1basename,'_PORCUPINE/',inputT1basename,'_ICD');
   % define filename prefixes
   rfp = strcat(rawFolder,'/',inputT1basename,'_');
   offp = strcat(outFolder,'/',inputT1basename,'_');
   oSTRAWffp = strcat(outSTRAWfolder,'/',inputT1basename,'_');
   oICDffp = strcat(outICDfolder,'/',inputT1basename,'_');
   % move final STRAW folder to output
   system(sprintf('mv -f %s %s',strawFolder,outSTRAWfolder));
   % move final ICD folder to output
   system(sprintf('mv -f %s %s',ntFolder,outICDfolder));
   % rename input T1 to RAW and move to output ICD folder
   system(sprintf('mv -f %s.nii %s.nii',strcat(rfp,'T1'),strcat(oICDffp,'RAW')));
   % move OVF to output ICD folder
   system(sprintf('mv -f %s.nii %s.nii',strcat(basedir,'/Processing/Stage1_StandardizeT1/',inputT1basename,'_OVF'),strcat(oICDffp,'OVF')));
   % copy SDN to T1 and move to output folder
   system(sprintf('cp -f %s.nii %s.nii',strcat(oSTRAWffp,'SDN'),strcat(offp,'T1')));
   % move visor targets file to output folder
   system(sprintf('mv -f %s.asc %s.asc',strcat(oICDffp,'VISOR_TARGETS'),strcat(offp,'VISOR_TARGETS')));
   status = true;
end