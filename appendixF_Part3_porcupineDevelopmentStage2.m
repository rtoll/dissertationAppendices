clear all
close all
clc

[~,userInfo] = system('id');
disp(userInfo)
[~,machineInfo] = system('hostname');
disp(machineInfo)
[~,bashRandom] = system('echo $RANDOM');
disp(bashRandom) % between 0 and 32767

pause(str2double(bashRandom)/32767*20);

addpath(genpath('/scratch/users/rtoll/matlab_functions_downloaded/'))
addpath(genpath('/scratch/users/rtoll/matlab_functions_custom/'))

baseDir = '/scratch/users/rtoll/PORCUPINE/';

inputT1s = dir(fullfile(baseDir,'Input'));
inputT1s = {inputT1s.name}';
inputT1s = inputT1s(contains(inputT1s,'.nii'));

% remove above in actual call stack

d = dir(fullfile(baseDir,'Output'));
d = {d.name}';
d = d(contains(d,'A') & ~contains(d,'.nii'));

while length(d) < 51

i=randperm(51,1);
volInputT1 = inputT1s{i};
basename = volInputT1(1:strfind(volInputT1,'.nii')-1);
subjOutputFolder = fullfile(baseDir,'Output',basename);
if ~exist(subjOutputFolder,'dir')
tic   
disp(i)
%
mkdir(subjOutputFolder)

%% phase02step01
% move inputT1 to output folder

% The input T1 is moved from the input folder to the subject's output
% folder. The system monitors the input folder to determine the processing
% queue with the while loop that monitors if the input folder is empty of
% files with a .nii or .nii.gz extension.
volRaw = fullfile(subjOutputFolder,strcat(basename,'_phase02step01_raw.nii'));
copyfile(fullfile(baseDir,'Input',volInputT1),volRaw) % change this to movefile

%% phase02step02
% convert NaNs (if present) to zero

% Not a number (NaN) values will affect downstream calculations and must be
% converted to zero.
volNoNans = fullfile(subjOutputFolder,strcat(basename,'_phase02step02_noNans.nii.gz'));
system(sprintf('fslmaths %s -nan %s',volRaw,volNoNans));

%% phase02step03
% reorient to standard (MNI) orientation (RL PA IS)

% The volume is then reoriented (not registered) to the standard
% MNI orientation (x = right to left, y = posterior to anterior,
% z = inferior to superior) by way of one or more 90 degree rotations.
volReorient = fullfile(subjOutputFolder,strcat(basename,'_phase02step03_reorient.nii.gz'));
system(sprintf('fslreorient2std %s %s',volNoNans,volReorient));

%% phase02step04
% radiological

% The volume is then verified to be in radiological (x = right to left)
% orientation and if not, the x axis is corrected to be so.
[~,currentOrientation] = system(sprintf('fslorient -getorient %s',volReorient));
volRadiological = fullfile(subjOutputFolder,strcat(basename,'_phase02step04_radiological.nii.gz'));
if strcmp(strtrim(currentOrientation),'NEUROLOGICAL')
   system(sprintf('fslswapdim %s -x y z %s',volReorient,volRadiological));
   system(sprintf('fslorient -forceradiological %s',volRadiological));
else
   copyfile(volReorient,volRadiological)
end

%% phase02step05
% make ovf

% The volume is then saved in its curent form as the original in visor
% format (OVF). Visor is neuronavigation software that presently requires
% the structural images it loads to be in INT16 format (using the fslhd
% command, one can inspect the data_type field). This is achieved by using 
% fslmaths to save the image with output data type (-odt) short. Saving the 
% OVF volume at this stage ensures a standardized version of the original volume
% is available for manual targeting should errors or failures occur in downstream
% processing.
volOVF = fullfile(subjOutputFolder,strcat(basename,'_phase02step05_ovf.nii.gz'));
system(sprintf('fslmaths %s %s -odt short',volRadiological,volOVF));

%% phase02step06
% set minimum value to 0 and convert to data type single (float)

% The volume's minimum value is then set to 0 to establish a common lower
% bound and the data type is then set to float. This is because the next
% step will involve interpolation and precision in calculation is
% significantly enhanced using float values instead of integer values.
[~,minVal] = system(sprintf('fslstats %s -p 0',volRadiological));
minVal = str2double(minVal);
volZeroFloat = fullfile(subjOutputFolder,strcat(basename,'_phase02step06_zeroFloat.nii.gz'));
system(sprintf('fslmaths %s -sub %u %s -odt float',volRadiological,minVal,volZeroFloat));

%% phase02step07
% make ISO

% The 
[~,xNumVox] = system(sprintf('fslval %s dim1',volZeroFloat));
xNumVox = str2double(xNumVox);
[~,xPixDim] = system(sprintf('fslval %s pixdim1',volZeroFloat));
xPixDim = str2double(xPixDim);
xMM = xNumVox * xPixDim;
xTrans = round((256 - xMM) / 2);
[~,yNumVox] = system(sprintf('fslval %s dim2',volZeroFloat));
yNumVox = str2double(yNumVox);
[~,yPixDim] = system(sprintf('fslval %s pixdim2',volZeroFloat));
yPixDim = str2double(yPixDim);
yMM = yNumVox * yPixDim;
yTrans = round((256 - yMM) / 2);
[~,zNumVox] = system(sprintf('fslval %s dim3',volZeroFloat));
zNumVox = str2double(zNumVox);
[~,zPixDim] = system(sprintf('fslval %s pixdim3',volZeroFloat));
zPixDim = str2double(zPixDim);
zMM = zNumVox * zPixDim;
zTrans = round((256 - zMM) / 2);
transToCenter = eye(4);
transToCenter(1:3,end) = [xTrans ; yTrans ; zTrans];
xform_zeroFloat_to_iso = fullfile(subjOutputFolder,strcat(basename,'_phase02step07a_xform_zeroFloat_to_iso.txt'));
dlmwrite(xform_zeroFloat_to_iso,transToCenter,' ');
volIso = fullfile(subjOutputFolder,strcat(basename,'_phase02step07b_iso.nii.gz'));
% make iso by flirting to template_zeros with xform_zeroFloat_to_iso
system(sprintf('flirt -in %s -ref %s/Admin/StructuralTemplates/templateZeros.nii.gz -applyxfm -init %s -out %s',volZeroFloat,baseDir,xform_zeroFloat_to_iso,volIso));

%% phase02step08
% robustfov

volGuillInitial = fullfile(subjOutputFolder,strcat(basename,'_phase02step08a_guillInitial.nii.gz'));
xform_guillInitial_to_guill = fullfile(subjOutputFolder,strcat(basename,'_phase02step08b_xform_guillInitial_to_guill.txt'));
system(sprintf('robustfov -i %s -m %s -r %s',volIso,xform_guillInitial_to_guill,volGuillInitial));
volGuill = fullfile(subjOutputFolder,strcat(basename,'_phase02step08c_guill.nii.gz'));
system(sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',volGuillInitial,volIso,volGuill,xform_guillInitial_to_guill));

%% phase02step09
% deep bet

volGuillBrain = fullfile(subjOutputFolder,strcat(basename,'_phase02step09a_guillBrain.nii.gz'));
system(sprintf('bet %s %s -f .75 -R',volGuill,volGuillBrain));
[~,guillBrainCoM] = system(sprintf('fslstats %s -C',volGuillBrain));
guillBrainCoM = str2num(guillBrainCoM);

volIsoBrainEyesCleanup = fullfile(subjOutputFolder,strcat(basename,'_phase02step09b_isoBrainEyesCleanup.nii.gz'));
system(sprintf('bet %s %s -m -f .1 -c %g %g %g -S',volIso,volIsoBrainEyesCleanup,guillBrainCoM(1),guillBrainCoM(2),guillBrainCoM(3)));
volIsoBrainEyesCleanupMask = fullfile(subjOutputFolder,strcat(basename,'_phase02step09c_isoBrainEyesCleanupMask.nii.gz'));
movefile(strrep(volIsoBrainEyesCleanup,'.nii.gz','_mask.nii.gz'),volIsoBrainEyesCleanupMask);
delete(strrep(volIsoBrainEyesCleanup,'.nii.gz','_skull.nii.gz'));

volIsoBrainNeckCleanup = fullfile(subjOutputFolder,strcat(basename,'_phase02step09d_isoBrainNeckCleanup.nii.gz'));
system(sprintf('bet %s %s -m -f .1 -c %g %g %g -B',volIso,volIsoBrainNeckCleanup,guillBrainCoM(1),guillBrainCoM(2),guillBrainCoM(3)));
volIsoBrainNeckCleanupMask = fullfile(subjOutputFolder,strcat(basename,'_phase02step09e_isoBrainNeckCleanupMask.nii.gz'));
movefile(strrep(volIsoBrainNeckCleanup,'.nii.gz','_mask.nii.gz'),volIsoBrainNeckCleanupMask);

volIsoBrainMask = fullfile(subjOutputFolder,strcat(basename,'_phase02step09f_isoBrainMask.nii.gz'));
system(sprintf('fslmaths %s -mas %s %s',volIsoBrainEyesCleanupMask,volIsoBrainNeckCleanupMask,volIsoBrainMask));
system(sprintf('fslmaths %s -s 3 %s -odt float',volIsoBrainMask,volIsoBrainMask));
system(sprintf('fslmaths %s -thr .5 -bin %s',volIsoBrainMask,volIsoBrainMask));
volIsoBrainMaskEroded = fullfile(subjOutputFolder,strcat(basename,'_phase02step09g_isoBrainMaskEroded.nii.gz'));
system(sprintf('fslmaths %s -ero -ero %s',volIsoBrainMask,volIsoBrainMaskEroded));
volIsoBrain = fullfile(subjOutputFolder,strcat(basename,'_phase02step09h_isoBrain.nii.gz'));
system(sprintf('fslmaths %s -mas %s %s -odt float',volIso,volIsoBrainMask,volIsoBrain));

%% phase02step10
% fast bias correction and simple segmentation

volIsoBrainBase = fullfile(subjOutputFolder,strcat(basename,'_phase02step10_isoBrain'));
system(sprintf('fast -I 10 -g --nopve -b -W 30 -O 10 -H .2 -o %s %s',volIsoBrainBase,volIsoBrain));
delete(strcat(volIsoBrainBase,'_seg.nii.gz'));
delete(strcat(volIsoBrainBase,'_seg_0.nii.gz'));
delete(strcat(volIsoBrainBase,'_seg_1.nii.gz'));
volIsoWhiteMatterMask = fullfile(subjOutputFolder,strcat(basename,'_phase02step10a_isoWhiteMatterMask.nii.gz'));
movefile(strcat(volIsoBrainBase,'_seg_2.nii.gz'),volIsoWhiteMatterMask);
volIsoWhiteMatterMaskClusters = fullfile(subjOutputFolder,strcat(basename,'_phase02step10b_isoWhiteMatterMaskClusters.nii.gz'));
system(sprintf('cluster -i %s -t 1 --connectivity=6 --no_table -o %s',volIsoWhiteMatterMask,volIsoWhiteMatterMaskClusters));
[~,maxClusterIndex] = system(sprintf('fslstats %s -p 100',volIsoWhiteMatterMaskClusters));
maxClusterIndex = str2double(maxClusterIndex);
volIsoWhiteMatterMaskMainCluster = fullfile(subjOutputFolder,strcat(basename,'_phase02step10c_isoWhiteMatterMaskMainCluster.nii.gz'));
system(sprintf('fslmaths %s -thr %u -bin %s',volIsoWhiteMatterMaskClusters,maxClusterIndex,volIsoWhiteMatterMaskMainCluster));
volIsoBrainBias = fullfile(subjOutputFolder,strcat(basename,'_phase02step10d_isoBrainBias.nii.gz'));
movefile(strcat(volIsoBrainBase,'_bias.nii.gz'),volIsoBrainBias);
volIsoBrainBiasMasked = fullfile(subjOutputFolder,strcat(basename,'_phase02step10e_isoBrainBiasMasked.nii.gz'));
system(sprintf('fslmaths %s -mas %s %s -odt float',volIsoBrainBias,volIsoBrainMaskEroded,volIsoBrainBiasMasked));
volIsoBrainBiasFill = fullfile(subjOutputFolder,strcat(basename,'_phase02step10f_isoBrainBiasFill.nii.gz'));
system(sprintf('fslsmoothfill -i %s -o %s -m %s',volIsoBrainBias,volIsoBrainBiasFill,volIsoBrainMask));
delete(strrep(volIsoBrainBiasFill,'.nii.gz','_idxmask.nii.gz'));
delete(strrep(volIsoBrainBiasFill,'.nii.gz','_init.nii.gz'));
delete(strrep(volIsoBrainBiasFill,'.nii.gz','_vol2.nii.gz'));
delete(strrep(volIsoBrainBiasFill,'.nii.gz','_vol32.nii.gz'));
volIsoBrainBiasSmooth = fullfile(subjOutputFolder,strcat(basename,'_phase02step10g_isoBrainBiasSmooth.nii.gz'));
system(sprintf('fslmaths %s -s 10 %s -odt float',volIsoBrainBiasFill,volIsoBrainBiasSmooth));
volIsoBrainBiasSmoothMasked = fullfile(subjOutputFolder,strcat(basename,'_phase02step10h_isoBrainBiasSmoothMasked.nii.gz'));
system(sprintf('fslmaths %s -binv -mul %s %s -odt float',volIsoBrainMaskEroded,volIsoBrainBiasSmooth,volIsoBrainBiasSmoothMasked));
volIsoBias = fullfile(subjOutputFolder,strcat(basename,'_phase02step10i_isoBias.nii.gz'));
system(sprintf('fslmaths %s -add %s %s -odt float',volIsoBrainBiasMasked,volIsoBrainBiasSmoothMasked,volIsoBias));
volIsoBiasCor = fullfile(subjOutputFolder,strcat(basename,'_phase02step10j_isoBiasCor.nii.gz'));
system(sprintf('fslmaths %s -div %s %s -odt float',volIso,volIsoBias,volIsoBiasCor));
volIsoBiasDiff = fullfile(subjOutputFolder,strcat(basename,'_phase02step10k_isoBiasDiff.nii.gz'));
system(sprintf('fslmaths %s -sub %s %s -odt float',volIsoBiasCor,volIso,volIsoBiasDiff));

%% phase02step11
% rescale

volIsoWhiteMatter = fullfile(subjOutputFolder,strcat(basename,'_phase02step11a_isoWhiteMatter.nii.gz'));
system(sprintf('fslmaths %s -mas %s %s -odt float',volIsoBiasCor,volIsoWhiteMatterMaskMainCluster,volIsoWhiteMatter));
[~,whiteMatterMedian] = system(sprintf('fslstats %s -P 50',volIsoWhiteMatter));
whiteMatterMedian = str2double(whiteMatterMedian);
threshMaxIntensity = 19999;
volRescale = fullfile(subjOutputFolder,strcat(basename,'_phase02step11b_rescale.nii.gz'));
[~,templateElecWhiteMatterMedian] = system(sprintf('fslstats %s/Admin/StructuralTemplates/templateElecWhiteMatter.nii.gz -P 50',baseDir));
templateElecWhiteMatterMedian = str2double(templateElecWhiteMatterMedian);
system(sprintf('fslmaths %s -div %g -mul %g %s -odt float',volIsoBiasCor,whiteMatterMedian,templateElecWhiteMatterMedian,volRescale));
system(sprintf('fslmaths %s -min %g %s -odt float',volRescale,threshMaxIntensity,volRescale));
volRescaleBrain = fullfile(subjOutputFolder,strcat(basename,'_phase02step11c_rescaleBrain.nii.gz'));
system(sprintf('fslmaths %s -mas %s %s -odt float',volRescale,volIsoBrainMask,volRescaleBrain));

%% phase02step12
% integer translation

volRescaleBrainRigid = fullfile(subjOutputFolder,strcat(basename,'_phase02step12a_rescaleBrainRigid.nii.gz'));
xform_rescaleBrain_to_rescaleBrainRigid = fullfile(subjOutputFolder,strcat(basename,'_phase02step12b_xform_rescaleBrain_to_rescaleBrainRigid.txt'));
system(sprintf('flirt -in %s -ref %s/Admin/StructuralTemplates/templateElecBrain.nii.gz -out %s -omat %s -cost normcorr -dof 6',volRescaleBrain,baseDir,volRescaleBrainRigid,xform_rescaleBrain_to_rescaleBrainRigid));

[~,rescaleBrainRigidCoM] = system(sprintf('fslstats %s -C',volRescaleBrainRigid));
rescaleBrainRigidCoM = str2num(rescaleBrainRigidCoM);
[~,rescaleBrainCoM] = system(sprintf('fslstats %s -C',volRescaleBrain));
rescaleBrainCoM = str2num(rescaleBrainCoM);
trans_rescaleBrain_to_translateBrain = eye(4);
trans_rescaleBrain_to_translateBrain(1:3,4) = round(rescaleBrainRigidCoM - rescaleBrainCoM)';

xform_rescaleBrain_to_translateBrain = fullfile(subjOutputFolder,strcat(basename,'_phase02step12c_xform_rescaleBrain_to_translateBrain.txt'));
dlmwrite(xform_rescaleBrain_to_translateBrain,trans_rescaleBrain_to_translateBrain,' ');

volTranslate = fullfile(subjOutputFolder,strcat(basename,'_phase02step12d_translate.nii.gz'));
system(sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',volRescale,volRescale,volTranslate,xform_rescaleBrain_to_translateBrain));
volTranslateBrain = fullfile(subjOutputFolder,strcat(basename,'_phase02step12e_translateBrain.nii.gz'));
system(sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',volRescaleBrain,volRescaleBrain,volTranslateBrain,xform_rescaleBrain_to_translateBrain));

%%

% rigidBrain
% cleanup brain
% flirtback

% rigid again for head
% applyxfm to head
% try 6/3mouth imgauss3sig
% flirt back

% double flirt to rigid last time
% mni dof 12
% do segs and fids
% apply elecs

% do measurements
% affine targets

% fnirt
% non lin targets

%% phase03step01
% rigid 1 and deep segmentation

volRigidBrainInitial = fullfile(subjOutputFolder,strcat(basename,'_phase03step01a_rigidBrainInitial.nii.gz'));
xform_translateBrain_to_rigidBrainInitial = fullfile(subjOutputFolder,strcat(basename,'_phase03step01b_xform_translateBrain_to_rigidBrainInitial.txt'));
system(sprintf('flirt -in %s -ref %s/Admin/StructuralTemplates/templateElecBrain.nii.gz -out %s -omat %s -cost normcorr -dof 6',volTranslateBrain,baseDir,volRigidBrainInitial,xform_translateBrain_to_rigidBrainInitial));
xform_rigidBrainInitial_to_translateBrain = fullfile(subjOutputFolder,strcat(basename,'_phase03step01c_xform_rigidBrainInitial_to_translateBrain.txt'));
system(sprintf('convert_xfm -omat %s -inverse %s',xform_rigidBrainInitial_to_translateBrain,xform_translateBrain_to_rigidBrainInitial));

volRigidBrainInitialBase = fullfile(subjOutputFolder,strcat(basename,'_phase03step01_rigidBrainInitial'));
system(sprintf('fast -N -W 30 -O 15 -H .2 -o %s %s',volRigidBrainInitialBase,volRigidBrainInitial));

delete(strcat(volRigidBrainInitialBase,'_mixeltype.nii.gz'));
delete(strcat(volRigidBrainInitialBase,'_pveseg.nii.gz'));
delete(strcat(volRigidBrainInitialBase,'_mixeltype.nii.gz'));

volRigidBrainInitialCSFpve = fullfile(subjOutputFolder,strcat(basename,'_phase03step01d_rigidBrainInitialCSFpve.nii.gz'));
movefile(strcat(volRigidBrainInitialBase,'_pve_0.nii.gz'),volRigidBrainInitialCSFpve);
volRigidBrainInitialGMpve = fullfile(subjOutputFolder,strcat(basename,'_phase03step01e_rigidBrainInitialGMpve.nii.gz'));
movefile(strcat(volRigidBrainInitialBase,'_pve_1.nii.gz'),volRigidBrainInitialGMpve);
volRigidBrainInitialWMpve = fullfile(subjOutputFolder,strcat(basename,'_phase03step01f_rigidBrainInitialWMpve.nii.gz'));
movefile(strcat(volRigidBrainInitialBase,'_pve_2.nii.gz'),volRigidBrainInitialWMpve);

nii = load_nii(volRigidBrainInitialCSFpve);
voxCSF = nii.img;
voxCSFhi = voxCSF >= 0;
voxCSFlo = voxCSF >= .5;

nii = load_nii(volRigidBrainInitialGMpve);
voxGM = nii.img;
voxGM = voxGM >= .5;
CC = bwconncomp(voxGM,6);
clusterLengths = cellfun(@length,CC.PixelIdxList);
voxGM = false(size(voxGM));
voxGM(CC.PixelIdxList{clusterLengths == max(clusterLengths)}) = true;

nii = load_nii(volRigidBrainInitialWMpve);
voxWM = nii.img;
voxWM = voxWM >= .5;
CC = bwconncomp(voxWM,6);
clusterLengths = cellfun(@length,CC.PixelIdxList);
voxWM = false(size(voxWM));
voxWM(CC.PixelIdxList{clusterLengths == max(clusterLengths)}) = true;

voxCM = voxGM | voxWM;

nii = load_nii(volRigidBrainInitial);
vox = nii.img;

voxGMvals = vox .* cast(voxGM,'like',vox);
voxGMvals = voxGMvals(voxGMvals > 0);
threshIntensity = mean(voxGMvals) - (2 * std(voxGMvals));
% thresh by mean(GM) - (2*SD(GM))

voxBin = vox >= threshIntensity;

% first cull in 3d
CC = bwconncomp(voxBin,6);
clusterLengths = cellfun(@length,CC.PixelIdxList);
voxBin = false(size(voxBin));
voxBin(CC.PixelIdxList{clusterLengths == max(clusterLengths)}) = true;

threshCrossSectionalArea = 5^2; % mm^2
for dim=1:3
   for d=1:256
      if dim == 1
         voxBinSlice = squeeze(voxBin(d,:,:));
      elseif dim == 2
         voxBinSlice = squeeze(voxBin(:,d,:));
      elseif dim == 3
         voxBinSlice = squeeze(voxBin(:,:,d));
      end
      CC = bwconncomp(voxBinSlice,4);
      voxBinSlice = false(size(voxBinSlice));
      for c=1:CC.NumObjects
         cluster = CC.PixelIdxList{c};
         if length(cluster) >= threshCrossSectionalArea
            voxBinSlice(cluster) = true;
         end      
      end
      if dim == 1
         voxBin(d,:,:) = voxBinSlice;
      elseif dim == 2
         voxBin(:,d,:) = voxBinSlice;
      elseif dim == 3
         voxBin(:,:,d) = voxBinSlice;
      end
   end
end

% second cull in 3d
CC = bwconncomp(voxBin,6);
clusterLengths = cellfun(@length,CC.PixelIdxList);
voxBin = false(size(voxBin));
voxBin(CC.PixelIdxList{clusterLengths == max(clusterLengths)}) = true;

% expand
voxBorder = false(size(voxBin));
voxBorder(1,:,:) = true;
voxBorder(end,:,:) = true;
voxBorder(:,1,:) = true;
voxBorder(:,end,:) = true;
voxBorder(:,:,1) = true;
voxBorder(:,:,end) = true;

for iterations=1:100
   disp(iterations) % tempdel
   startSum = sum(voxBin(:));
   voxNull = ~voxBin;
   voxAdd = bwmorph3(voxNull,'remove') & ~voxBorder;
   voxBin = voxBin | (voxAdd & ~voxCSFhi & voxCM);
   if sum(voxBin(:)) == startSum
      break
   end
end
for iterations=1:3
   voxNull = ~voxBin;
   voxAdd = bwmorph3(voxNull,'remove') & ~voxBorder;
   voxBin = voxBin | (voxAdd & ~voxCSFlo & voxCM);
   voxBin = voxBin & bwmorph3(voxBin,'majority');
end
for iterations=1:2
   voxNull = ~voxBin;
   voxAdd = bwmorph3(voxNull,'remove') & ~voxBorder;
   voxBin = voxBin | voxAdd;
end
voxBin = bwmorph3(voxBin,'majority');
voxBin = imfill(voxBin,6,'holes');
volRigidBrainInitialRefined = fullfile(subjOutputFolder,strcat(basename,'_phase03step01g_rigidBrainInitialRefined.nii.gz'));
nii.img = vox .* cast(voxBin,'like',vox);
save_nii(nii,volRigidBrainInitialRefined);

volTranslateBrainRefined = fullfile(subjOutputFolder,strcat(basename,'_phase03step01h_translateBrainRefined.nii.gz'));
system(sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',volRigidBrainInitialRefined,volRigidBrainInitialRefined,volTranslateBrainRefined,xform_rigidBrainInitial_to_translateBrain));

%% phase03step02
% native xform and head cleaning

volAffineBrain = fullfile(subjOutputFolder,strcat(basename,'_phase03step02a_affineBrain.nii.gz'));
xform_translateBrainRefined_to_affineBrain = fullfile(subjOutputFolder,strcat(basename,'_phase03step02b_xform_translateBrainRefined_to_affineBrain.txt'));
system(sprintf('flirt -in %s -ref %s/Admin/StructuralTemplates/templateElecBrain.nii.gz -out %s -omat %s -cost normcorr',volTranslateBrainRefined,baseDir,volAffineBrain,xform_translateBrainRefined_to_affineBrain));

volNativeBrain = fullfile(subjOutputFolder,strcat(basename,'_phase03step02c_nativeBrain.nii.gz'));
xform_translateBrainRefined_to_nativeBrain = fullfile(subjOutputFolder,strcat(basename,'_phase03step02d_xform_translateBrainRefined_to_nativeBrain.txt'));
system(sprintf('flirt -in %s -ref %s -out %s -omat %s -cost normcorr -dof 6',volTranslateBrainRefined,volAffineBrain,volNativeBrain,xform_translateBrainRefined_to_nativeBrain));
xform_nativeBrain_to_translateBrainRefined = fullfile(subjOutputFolder,strcat(basename,'_phase03step02e_xform_nativeBrain_to_translateBrainRefined.txt'));
system(sprintf('convert_xfm -omat %s -inverse %s',xform_nativeBrain_to_translateBrainRefined,xform_translateBrainRefined_to_nativeBrain));

volNativePreclean = fullfile(subjOutputFolder,strcat(basename,'_phase03step02f_nativePreclean.nii.gz'));
system(sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',volTranslate,volTranslate,volNativePreclean,xform_translateBrainRefined_to_nativeBrain));

[~,nativeBrainCoM] = system(sprintf('fslstats %s -C',volNativeBrain));
nativeBrainCoM = str2num(nativeBrainCoM);
volNativeBrainBase = fullfile(subjOutputFolder,strcat(basename,'_phase03step02_nativeBrain'));
system(sprintf('bet %s %s -c %g %g %g -A',volNativePreclean,volNativeBrainBase,nativeBrainCoM(1),nativeBrainCoM(2),nativeBrainCoM(3)));

volNativeOutskullMask = fullfile(subjOutputFolder,strcat(basename,'_phase03step02g_nativeOutskullMask.nii.gz'));
movefile(strcat(volNativeBrainBase,'_outskull_mask.nii.gz'),volNativeOutskullMask);

volNativeOutskinMesh = fullfile(subjOutputFolder,strcat(basename,'_phase03step02h_nativeOutskinMesh.nii.gz'));
movefile(strcat(volNativeBrainBase,'_outskin_mesh.nii.gz'),volNativeOutskinMesh);

delete(strcat(volNativeBrainBase,'_outskin_mask.nii.gz'));
delete(strcat(volNativeBrainBase,'_outskin_mesh.vtk'));
delete(strcat(volNativeBrainBase,'_outskull_mesh.nii.gz'));
delete(strcat(volNativeBrainBase,'_outskull_mesh.vtk'));
delete(strcat(volNativeBrainBase,'_skull_mask.nii.gz'));
delete(strcat(volNativeBrainBase,'.nii.gz'));
delete(strcat(volNativeBrainBase,'_inskull_mask.nii.gz'));
delete(strcat(volNativeBrainBase,'_inskull_mesh.nii.gz'));
delete(strcat(volNativeBrainBase,'_inskull_mesh.vtk'));
delete(strcat(volNativeBrainBase,'_mesh.vtk'));







nii = load_nii(volNativeOutskullMask);
voxNativeOutskullMask = logical(nii.img);

nii = load_nii(volNativeOutskinMesh);
voxNativeOutskinMesh = logical(nii.img);

voxNativeOutskullMaskDilate = voxNativeOutskullMask;
while isempty(find(voxNativeOutskullMaskDilate & voxNativeOutskinMesh,1))
   voxNativeOutskullMaskDilate = voxNativeOutskullMaskDilate | bwmorph3(~voxNativeOutskullMaskDilate,'remove') & ~voxBorder;
end
nii.img = voxNativeOutskullMaskDilate;
volNativeOutskullMaskDilated = fullfile(subjOutputFolder,strcat(basename,'_phase03step02i_nativeOutskullMask.nii.gz'));
save_nii(nii,volNativeOutskullMaskDilated);

volNativeGauss = fullfile(subjOutputFolder,strcat(basename,'_phase03step02j_nativeGauss.nii.gz'));
system(sprintf('fslmaths %s -s 3 %s -odt float',volNativePreclean,volNativeGauss));

volNativeGaussSubtracted = fullfile(subjOutputFolder,strcat(basename,'_phase03step02k_nativeGaussSubtracted.nii.gz'));
system(sprintf('fslmaths %s -sub %s %s -odt float',volNativePreclean,volNativeGauss,volNativeGaussSubtracted));

voxNativeGaussSubtracted = load_nii(volNativeGaussSubtracted);
voxNativePreclean = load_nii(volNativePreclean);
threshMinSubtractedIntensity = 250;
threshMinPrecleanIntensity = 500;
voxNativeGaussSubtracted = voxNativeGaussSubtracted >= threshMinSubtractedIntensity & voxNativePreclean >= threshMinPrecleanIntensity;
voxNativeGaussSubtracted = voxNativeGaussSubtracted | voxNativeOutskullMaskDilate;

% now ready to cull

nii = load_nii(volNativeBrain);
voxNativeBrain = nii.img;
[brainX,brainY,brainZ] = ind2sub(size(voxNativeBrain),find(voxNativeBrain > 0));












%%
system(sprintf('chmod -R 777 %s',subjOutputFolder));
toc/60


end

d = dir(fullfile(baseDir,'Output'));
d = {d.name}';
d = d(contains(d,'A') & ~contains(d,'.nii'));
end

%%

% nii = load_nii(volNativeBrain);
% voxBrain = nii.img;
% [brainX,brainY,brainZ] = ind2sub(size(voxBrain),find(voxBrain > 0));



% nii = load_nii(volNativePreclean);
% vox = nii.img;
% threshMinIntensity = 500;
% threshOffset = 2; % voxels
% vox(vox < threshMinIntensity) = threshMinIntensity;
% voxContrast26 = single(zeros(256,256,256,26));
% transKey = NaN(26,3);
% i=1;
% for xx=-1:1
%    for yy=-1:1
%       for zz=-1:1
%          trans = [xx,yy,zz];
%          if sum(abs(trans)) ~= 0
%             if sum(abs(trans)) == 1
%                % face
%                transAdj = round(trans .* threshOffset ./ sqrt(1));
%             elseif sum(abs(trans)) == 2
%                transAdj = round(trans .* threshOffset ./ sqrt(2));
%                % edge
%             elseif sum(abs(trans)) == 3
%                transAdj = round(trans .* threshOffset ./ sqrt(3));
%                % point
%             end
%             voxContrast26(:,:,:,i) = vox ./ voxNN(vox,[transAdj(1),transAdj(2),transAdj(3)],threshMinIntensity);
%             transKey(i,:) = trans;
%             i=i+1;
%          end
%       end
%    end
% end
% voxContrastOctants = single(zeros(size(vox)));
% voxLeft = false(size(vox));
% voxLeft(1:min(brainX) + floor(.5 * range(brainX)),:,:) = true;
% voxPosterior = false(size(vox));
% voxPosterior(:,1:min(brainY) + floor(.5 * range(brainY)),:) = true;
% voxInferior = false(size(vox));
% voxInferior(:,:,1:min(brainZ)-1) = true;
% for octant=1:8
%    if octant == 1
%       % LPI
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) <= 0 & transKey(:,2) <= 0 & transKey(:,3) <= 0),[],4);
%       voxContrastOctants(voxLeft & voxPosterior & voxInferior) = voxTemp(voxLeft & voxPosterior & voxInferior);
%    elseif octant == 2
%       % RPI
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) >= 0 & transKey(:,2) <= 0 & transKey(:,3) <= 0),[],4);
%       voxContrastOctants(~voxLeft & voxPosterior & voxInferior) = voxTemp(~voxLeft & voxPosterior & voxInferior);
%    elseif octant == 3
%       % LAI
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) <= 0 & transKey(:,2) >= 0 & transKey(:,3) <= 0),[],4);
%       voxContrastOctants(voxLeft & ~voxPosterior & voxInferior) = voxTemp(voxLeft & ~voxPosterior & voxInferior);
%    elseif octant == 4
%       % RAI
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) >= 0 & transKey(:,2) >= 0 & transKey(:,3) <= 0),[],4);
%       voxContrastOctants(~voxLeft & ~voxPosterior & voxInferior) = voxTemp(~voxLeft & ~voxPosterior & voxInferior);
%    elseif octant == 5
%       % LPS
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) <= 0 & transKey(:,2) <= 0 & transKey(:,3) >= 0),[],4);
%       voxContrastOctants(voxLeft & voxPosterior & ~voxInferior) = voxTemp(voxLeft & voxPosterior & ~voxInferior);
%    elseif octant == 6
%       % RPS
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) >= 0 & transKey(:,2) <= 0 & transKey(:,3) >= 0),[],4);
%       voxContrastOctants(~voxLeft & voxPosterior & ~voxInferior) = voxTemp(~voxLeft & voxPosterior & ~voxInferior);
%    elseif octant == 7
%       % LAS
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) <= 0 & transKey(:,2) >= 0 & transKey(:,3) >= 0),[],4);
%       voxContrastOctants(voxLeft & ~voxPosterior & ~voxInferior) = voxTemp(voxLeft & ~voxPosterior & ~voxInferior);
%    elseif octant == 8
%       % RAS
%       voxTemp = max(voxContrast26(:,:,:,transKey(:,1) >= 0 & transKey(:,2) >= 0 & transKey(:,3) >= 0),[],4);
%       voxContrastOctants(~voxLeft & ~voxPosterior & ~voxInferior) = voxTemp(~voxLeft & ~voxPosterior & ~voxInferior);
%    end
% end
% nii.img = voxContrastOctants;
% volNativeContrastOctants = fullfile(subjOutputFolder,strcat(basename,'_phase03step02j_nativeContrastOctants.nii.gz'));
% save_nii(nii,volNativeContrastOctants);
