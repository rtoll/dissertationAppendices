function [status] = porcupineApplyElectrodes(inputT1basename,basedir,fsldir)
   status = false;
   setenv('FSLDIR',fsldir);
   setenv('FSLOUTPUTTYPE', 'NIFTI');
   % define and create directories
   strawFolder = strcat(basedir,'/Processing/Stage2_StructuralRegistrationAndWarps/',inputT1basename,'_STRAW');
   ntFolder = strcat(basedir,'/Processing/Stage3_ApplyTemplates/',inputT1basename,'_TEMPLATES'); % native templates folder
   % define filenames
   sffp = strcat(strawFolder,'/',inputT1basename,'_'); % strawFolder with filename prefix for strcat
   sdn = strcat(sffp,'SDN'); % standardized and denoised
   ntffp = strcat(ntFolder,'/',inputT1basename,'_'); % ntFolder with filename prefix for strcat
   % load refined coordinates for easyCap and brainCap
   load(strcat(basedir,'/Admin/easyCap.mat')) % cohen cap
   load(strcat(basedir,'/Admin/waveguardSubsetCap.mat')) % older cap
   % load master target labels and IDs
   load(strcat(basedir,'/Admin/masterTargets.mat'))
   % obtain T1 dimensions
   [~,pixdim1] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',sdn,'pixdim1'));
   pixdim1 = str2double(pixdim1);
   [~,pixdim2] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',sdn,'pixdim2'));
   pixdim2 = str2double(pixdim2);
   [~,pixdim3] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',sdn,'pixdim3'));
   pixdim3 = str2double(pixdim3);
   voxToMM = [pixdim1,pixdim2,pixdim3];
   % define threshold
   threshTapeDist = 12; % mm, the interval at which the virtual tape "bends" while measuring head
   % load sdn
   t = porcupineNiiToVoxStan(sdn,fsldir);
   tb = logical(t);
   load(strcat(ntffp,'VTM.mat')) % VTM is in neuro voxel coordinates LRPAISvx
   % define fiducials
   nudgeElecMM = 1; % offset from head
   nudgeElecVX = round(nudgeElecMM / sqrt(sum(voxToMM .^ 2)));
   nudgeTapeMM = 1; % offset from head
   nudgeTapeVX = round(nudgeTapeMM / sqrt(sum(voxToMM .^ 2)));
   earSearchMM = 30;
   earSearchVX = round(earSearchMM / pixdim2);
   NASvx = vtm(vtm(:,1) == 210,2:end);
   INIvx = vtm(vtm(:,1) == 260,2:end);
   LPAvx = vtm(vtm(:,1) == 221,2:end);
   RPAvx = vtm(vtm(:,1) == 222,2:end);
   MPMvx = (NASvx + INIvx) ./ 2;
   MPEvx = (LPAvx + RPAvx) ./ 2;
   ORGvx = (MPMvx + MPEvx) ./ 2;
   CPvx = cross( (RPAvx - LPAvx) , (NASvx - INIvx) );
   CPvx = CPvx(end,:);
   CZvx = porcupineBresenhamFirstContactPlusNudge(tb,ORGvx,CPvx,1,[NASvx;INIvx;LPAvx;RPAvx]);
   FIDSvx = [NASvx;INIvx;LPAvx;RPAvx;CZvx;ORGvx;MPMvx;MPEvx];
   FIDSmm = FIDSvx .* repmat(voxToMM,[8 1]);
   stretchX = porcupineDistbwpts(FIDSmm(3,:),FIDSmm(4,:)) / 2;
   stretchY = porcupineDistbwpts(FIDSmm(1,:),FIDSmm(2,:)) / 2;
   stretchZ = porcupineDistbwpts(FIDSmm(5,:),FIDSmm(6,:));
   stretch = [stretchX,stretchY,stretchZ];
   %% write initial visor standard targets asc in case of quarantine
   abridgedTargets = sortrows(masterTargets,4);
   j=1;
   for i=1:size(abridgedTargets,1)
      if ~isnan(abridgedTargets{i,4}) && abridgedTargets{i,1} < 1000
         abridgedTargetsName{j,1} = abridgedTargets{i,3};
         abridgedTargetsVX(j,1:3) = vtm(vtm(:,1) == abridgedTargets{i,1},2:end);
         j=j+1;
      end
   end
   abridgedTargetsNE = porcupineCoordXfmNeuroVoxToVisorNE(NASvx,LPAvx,RPAvx,[pixdim1,pixdim2,pixdim3],[size(t,1),size(t,2),size(t,3)],abridgedTargetsVX);
   preamble{1} = sprintf('#\n');
   preamble{2} = sprintf('# use the following MRI coords for fiducials:\n');
   preamble{3} = sprintf('# Nasion = [%u %u %u]\n',NASvx(2),size(t,1) - NASvx(1) + 1,NASvx(3));
   preamble{4} = sprintf('# Left Ear = [%u %u %u]\n',LPAvx(2),size(t,1) - LPAvx(1) + 1,LPAvx(3));
   preamble{5} = sprintf('# Right Ear = [%u %u %u]\n',RPAvx(2),size(t,1) - RPAvx(1) + 1,RPAvx(3));
   preamble{6} = sprintf('# \n');
   preamble{7} = sprintf('# \n');
   preamble{8} = sprintf('# porcupine generated targets, direct questions/comments to russ, 512-663-3142, rtoll@stanford.edu\n');
   preamble{9} = sprintf('#\n');
   preamble{10} = sprintf('# %u position markers\n',length(abridgedTargetsName));
   preamble{11} = sprintf('NumberPositions=	%u\n',length(abridgedTargetsName));
   preamble{12} = sprintf('UnitPosition	mm\n');
   preamble{13} = sprintf('HSPTransformed	false\n');
   preamble{14} = sprintf('Positions\n');
   fclose('all');
   fid = fopen(strcat(ntffp,'VISOR_TARGETS.asc'),'w');
   for i=1:size(preamble,2)
      fwrite(fid,preamble{i});
   end
   labels = '';
   for i=1:size(abridgedTargetsNE,1)
      fwrite(fid,sprintf('%s :\t%g\t%g\t%g\n',abridgedTargetsName{i},abridgedTargetsNE(i,1),abridgedTargetsNE(i,2),abridgedTargetsNE(i,3)));
      labels = sprintf('%s\t%s',labels,abridgedTargetsName{i});
   end
   fwrite(fid,sprintf('Labels\n'));
   fwrite(fid,sprintf('%s\n',labels(2:end)));
   fclose('all');
   clear abridged* preamble labels
   %% bracketing
   fidsCap = [0,1,0;0,-1,0;-1,0,0;1,0,0];
   cap = fidsCap;
   capXfm = cap .* repmat(stretch,[length(cap) 1]) + repmat(FIDSmm(6,:),[length(cap) 1]);
   for k=0:5
      if k == 0
         angLo = -45;
         angHi = 45;
         ang = zeros(91,1);
         val = zeros(91,1);
      else
         ang = zeros(21,1);
         val = zeros(21,1);
      end
      j=1;
      for alpha = angLo : 10^-k : angHi
         capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(6,:),alpha,0,0);
         nasZdist = capXfmR(1,3) - FIDSmm(1,3);
         iniZdist = capXfmR(2,3) - FIDSmm(2,3);
         ang(j) = alpha;
         val(j) = abs(nasZdist - iniZdist);
         j=j+1;
      end
      [~,troughIndices] = findpeaks(-val);
      troughIndicesDistFromCenter = troughIndices - (j / 2);
      troughIndex = troughIndices(find(troughIndicesDistFromCenter == min(troughIndicesDistFromCenter),1,'first'));
      if isempty(troughIndex)
         troughIndex = j / 2;
      end
      angLo = ang(troughIndex - 1);
      angMin = ang(troughIndex);
      angHi = ang(troughIndex + 1);
   end
   alpha = angMin;
   for k=0:5
      if k == 0
         angLo = -45;
         angHi = 45;
         ang = zeros(91,1);
         val = zeros(91,1);
      else
         ang = zeros(21,1);
         val = zeros(21,1);
      end
      j=1;
      for beta = angLo : 10^-k : angHi
         capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(6,:),alpha,beta,0);
         lpaZdist = capXfmR(3,3) - FIDSmm(3,3);
         rpaZdist = capXfmR(4,3) - FIDSmm(4,3);
         ang(j) = beta;
         val(j) = abs(lpaZdist - rpaZdist);
         j=j+1;
      end
      [~,troughIndices] = findpeaks(-val);
      troughIndicesDistFromCenter = troughIndices - (j / 2);
      troughIndex = troughIndices(find(troughIndicesDistFromCenter == min(troughIndicesDistFromCenter),1,'first'));
      if isempty(troughIndex)
         troughIndex = j / 2;
      end
      angLo = ang(troughIndex - 1);
      angMin = ang(troughIndex);
      angHi = ang(troughIndex + 1);
   end
   beta = angMin;
   for k=0:5
      if k == 0
         angLo = -45;
         angHi = 45;
         ang = zeros(91,1);
         val = zeros(91,1);
      else
         ang = zeros(21,1);
         val = zeros(21,1);
      end
      j=1;
      for gamma = angLo : 10^-k : angHi
         capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(6,:),alpha,beta,gamma);
         lpaYdist = capXfmR(3,2) - FIDSmm(3,2);
         rpaYdist = capXfmR(4,2) - FIDSmm(4,2);
         ang(j) = gamma;
         val(j) = abs(lpaYdist - rpaYdist);
         j=j+1;
      end
      [~,troughIndices] = findpeaks(-val);
      troughIndicesDistFromCenter = troughIndices - (j / 2);
      troughIndex = troughIndices(find(troughIndicesDistFromCenter == min(troughIndicesDistFromCenter),1,'first'));
      if isempty(troughIndex)
         troughIndex = j / 2;
      end
      angLo = ang(troughIndex - 1);
      angMin = ang(troughIndex);
      angHi = ang(troughIndex + 1);
   end
   gammaEquator = angMin;
   for k=0:5
      if k == 0
         angLo = -45;
         angHi = 45;
         ang = zeros(91,1);
         val = zeros(91,1);
      else
         ang = zeros(21,1);
         val = zeros(21,1);
      end
      j=1;
      for gamma = angLo : 10^-k : angHi
         capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(6,:),alpha,beta,gamma);
         nasXdist = capXfmR(1,1) - FIDSmm(1,1);
         iniXdist = capXfmR(2,1) - FIDSmm(2,1);
         ang(j) = gamma;
         val(j) = abs(nasXdist - iniXdist);
         j=j+1;
      end
      [~,troughIndices] = findpeaks(-val);
      troughIndicesDistFromCenter = troughIndices - (j / 2);
      troughIndex = troughIndices(find(troughIndicesDistFromCenter == min(troughIndicesDistFromCenter),1,'first'));
      if isempty(troughIndex)
         troughIndex = j / 2;
      end
      angLo = ang(troughIndex - 1);
      angMin = ang(troughIndex);
      angHi = ang(troughIndex + 1);
   end
   gammaMeridian = angMin;
   gamma = mean([gammaEquator,gammaMeridian]);
   rotAngs = [alpha,beta,gamma,gammaMeridian,gammaEquator];
   save(strcat(ntffp,'ROTANGS.mat'),'rotAngs','FIDSvx','FIDSmm','voxToMM','stretch')
   %% measure head
   theta = linspace(0,180,1000);
   meridianCap(:,1) = zeros(1,1000);
   meridianCap(:,2) = cosd(theta);
   meridianCap(:,3) = sind(theta);
   equatorCap(:,1) = cosd(theta);
   equatorCap(:,2) = zeros(1,1000);
   equatorCap(:,3) = sind(theta);
   theta = linspace(0,360,1000);
   circumferenceCap(:,1) = cosd(theta);
   circumferenceCap(:,2) = sind(theta);
   circumferenceCap(:,3) = tand(20) .* ones(1,1000);
   clear cap*
   cap = meridianCap;
   capXfm = cap .* repmat(stretch,[length(cap) 1]) + repmat(FIDSmm(7,:),[length(cap) 1]);
   capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(7,:),alpha,beta,gammaMeridian);
   capXfmRE = 2 * (capXfmR - repmat(FIDSmm(7,:),[length(cap) 1])) + repmat(FIDSmm(7,:),[length(cap) 1]);
   capXfmREvx = capXfmRE ./ repmat(voxToMM,[length(cap) 1]);
   for i=1:length(cap)
      capFCvx(i,:) = porcupineBresenhamFirstContactPlusNudge(tb,MPMvx,capXfmREvx(i,:),nudgeTapeVX + 2,FIDSvx);
   end
   meridianCapFCvx = capFCvx;
   capFCmm = capFCvx .* repmat(voxToMM,[size(capFCvx,1),1]);
   index = 1;
   j=1;
   while ~isempty(index)
      startPoint = capFCmm(index,:);
      j=j+1;
      d = zeros(size(capFCmm,1),1);
      for i=index:size(capFCmm,1)
         d(i,1) = porcupineDistbwpts(startPoint,capFCmm(i,:));
      end
      index = find(d > threshTapeDist,1,'first');
   end
   NAStoINIcm = ((j-1) * threshTapeDist + round(porcupineDistbwpts(startPoint,capFCmm(end,:)))) / 10;
   clear cap*
   cap = equatorCap;
   capXfm = cap .* repmat(stretch,[length(cap) 1]) + repmat(FIDSmm(8,:),[length(cap) 1]);
   capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(7,:),alpha,beta,gammaEquator);
   capXfmRE = 2 * (capXfmR - repmat(FIDSmm(8,:),[length(cap) 1])) + repmat(FIDSmm(8,:),[length(cap) 1]);
   capXfmREvx = capXfmRE ./ repmat(voxToMM,[length(cap) 1]);
   for i=1:length(cap)
      capFCvx(i,:) = porcupineBresenhamFirstContactPlusNudge(tb,MPEvx,capXfmREvx(i,:),nudgeTapeVX + 1,FIDSvx);
   end
   equatorCapFCvx = capFCvx;
   capFCmm = capFCvx .* repmat(voxToMM,[size(capFCvx,1),1]);
   index = 1;
   j=1;
   while ~isempty(index)
      startPoint = capFCmm(index,:);
      j=j+1;
      d = zeros(size(capFCmm,1),1);
      for i=index:size(capFCmm,1)
         d(i,1) = porcupineDistbwpts(startPoint,capFCmm(i,:));
      end
      index = find(d > threshTapeDist,1,'first');
   end
   LPAtoRPAcm = ((j-1) * threshTapeDist + round(porcupineDistbwpts(startPoint,capFCmm(end,:)))) / 10;   
   clear cap*
   cap = circumferenceCap;
   capXfm = cap .* repmat(stretch,[length(cap) 1]) + repmat(FIDSmm(6,:),[length(cap) 1]);
   capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(6,:),alpha,beta,gamma);
   capXfmRE = 2 * (capXfmR - repmat(FIDSmm(6,:),[length(cap) 1])) + repmat(FIDSmm(6,:),[length(cap) 1]);
   capXfmREvx = capXfmRE ./ repmat(voxToMM,[length(cap) 1]);
   for i=1:length(cap)
      capFCvx(i,:) = porcupineBresenhamFirstContactPlusNudge(tb,ORGvx,capXfmREvx(i,:),nudgeTapeVX,FIDSvx);
   end
   circumferenceCapFCvx = capFCvx;
   capFCmm = capFCvx .* repmat(voxToMM,[size(capFCvx,1),1]);
   index = 1;
   j=1;
   while ~isempty(index)
      startPoint = capFCmm(index,:);
      j=j+1;
      d = zeros(size(capFCmm,1),1);
      for i=index:size(capFCmm,1)
         d(i,1) = porcupineDistbwpts(startPoint,capFCmm(i,:));
      end
      index = find(d > threshTapeDist,1,'first');
   end
   CIRCcm = ((j-1) * threshTapeDist + round(porcupineDistbwpts(startPoint,capFCmm(end,:)))) / 10;
   %% apply transforms 
   clear cap*
   cap = easyCap;
   capXfm = cap .* repmat(stretch,[length(cap) 1]) + repmat(FIDSmm(6,:),[length(cap) 1]);
   capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(6,:),alpha,beta,gamma);
   capXfmRE = 2 * (capXfmR - repmat(FIDSmm(6,:),[length(cap) 1])) + repmat(FIDSmm(6,:),[length(cap) 1]);
   capXfmREvx = capXfmRE ./ repmat(voxToMM,[length(cap) 1]);
   for i=1:length(cap)
      capFCvx(i,:) = porcupineBresenhamFirstContactPlusNudge(tb,ORGvx,capXfmREvx(i,:),nudgeElecVX,FIDSvx);
   end
   % check ears, easyCap = 53,47
   j=1;
   for i=0:earSearchVX
      capLEvxEndPoint = [capXfmREvx(53,1),capXfmREvx(53,2) - i,capXfmREvx(53,3)];
      capLEvx(j,:) = porcupineBresenhamFirstContactPlusNudge(tb,ORGvx,capLEvxEndPoint,nudgeElecVX,FIDSvx);
      capLEdist(j,:) = porcupineDistbwpts(ORGvx,capLEvx(j,:));
      j=j+1;
   end
   capFCvx(53,:) = capLEvx(find(capLEdist == min(capLEdist),1,'first'),:);
   j=1;
   for i=0:earSearchVX
      capREvxEndPoint = [capXfmREvx(47,1),capXfmREvx(47,2) - i,capXfmREvx(47,3)];
      capREvx(j,:) = porcupineBresenhamFirstContactPlusNudge(tb,ORGvx,capREvxEndPoint,nudgeElecVX,FIDSvx);
      capREdist(j,:) = porcupineDistbwpts(ORGvx,capREvx(j,:));
      j=j+1;
   end
   capFCvx(47,:) = capREvx(find(capREdist == min(capREdist),1,'first'),:);
   capFCmm = capFCvx .* repmat(voxToMM,[length(cap) 1]);
   easyCapFCvx = capFCvx;
   easyCapFCmm = capFCmm;
   clear cap*
   cap = waveguardSubsetCap;
   capXfm = cap .* repmat(stretch,[length(cap) 1]) + repmat(FIDSmm(6,:),[length(cap) 1]);
   capXfmR = porcupineSequentialRotateAboutAxis(capXfm,1,FIDSmm(6,:),alpha,beta,gamma);
   capXfmRE = 2 * (capXfmR - repmat(FIDSmm(6,:),[length(cap) 1])) + repmat(FIDSmm(6,:),[length(cap) 1]);
   capXfmREvx = capXfmRE ./ repmat(voxToMM,[length(cap) 1]);
   for i=1:length(cap)
      capFCvx(i,:) = porcupineBresenhamFirstContactPlusNudge(tb,ORGvx,capXfmREvx(i,:),nudgeElecVX,FIDSvx);
   end
   capFCmm = capFCvx .* repmat(voxToMM,[length(cap) 1]);
   waveguardSubsetCapFCvx = capFCvx;
   % append electrode positions to vtm and overwrite
   % 1000 series is easyCap
   for i=1:length(easyCapFCvx)
      vtm(end+1,1) = 1000 + i;
      vtm(end,2) = easyCapFCvx(i,1);
      vtm(end,3) = easyCapFCvx(i,2);
      vtm(end,4) = easyCapFCvx(i,3);
   end
   % 1100 series is waveguardSubset
   for i=1:length(waveguardSubsetCapFCvx)
      vtm(end+1,1) = 1100 + i;
      vtm(end,2) = waveguardSubsetCapFCvx(i,1);
      vtm(end,3) = waveguardSubsetCapFCvx(i,2);
      vtm(end,4) = waveguardSubsetCapFCvx(i,3);
   end
   save(strcat(ntffp,'VTM.mat'),'vtm')
   %% write visor standard targets asc
   abridgedTargets = sortrows(masterTargets,4);
   j=1;
   for i=1:size(abridgedTargets,1)
      if ~isnan(abridgedTargets{i,4})
         abridgedTargetsName{j,1} = abridgedTargets{i,3};
         abridgedTargetsVX(j,1:3) = vtm(vtm(:,1) == abridgedTargets{i,1},2:end);
         j=j+1;
      end
   end
   abridgedTargetsNE = porcupineCoordXfmNeuroVoxToVisorNE(NASvx,LPAvx,RPAvx,[pixdim1,pixdim2,pixdim3],[size(t,1),size(t,2),size(t,3)],abridgedTargetsVX);
   preamble{1} = sprintf('#\n');
   preamble{2} = sprintf('# use the following MRI coords for fiducials:\n');
   preamble{3} = sprintf('# Nasion = [%u %u %u]\n',NASvx(2),size(t,1) - NASvx(1) + 1,NASvx(3));
   preamble{4} = sprintf('# Left Ear = [%u %u %u]\n',LPAvx(2),size(t,1) - LPAvx(1) + 1,LPAvx(3));
   preamble{5} = sprintf('# Right Ear = [%u %u %u]\n',RPAvx(2),size(t,1) - RPAvx(1) + 1,RPAvx(3));
   preamble{6} = sprintf('# estimated head measurements:\n');
   preamble{7} = sprintf('# circ = %g cm, nas to ini = %g cm, tragus to tragus = %g cm\n',CIRCcm,NAStoINIcm,LPAtoRPAcm);
   preamble{8} = sprintf('# porcupine generated targets, direct questions/comments to russ, 512-663-3142, rtoll@stanford.edu\n');
   preamble{9} = sprintf('#\n');
   preamble{10} = sprintf('# %u position markers\n',length(abridgedTargetsName));
   preamble{11} = sprintf('NumberPositions=	%u\n',length(abridgedTargetsName));
   preamble{12} = sprintf('UnitPosition	mm\n');
   preamble{13} = sprintf('HSPTransformed	false\n');
   preamble{14} = sprintf('Positions\n');
   fclose('all');
   fid = fopen(strcat(ntffp,'VISOR_TARGETS.asc'),'w');
   for i=1:size(preamble,2)
      fwrite(fid,preamble{i});
   end
   labels = '';
   for i=1:size(abridgedTargetsNE,1)
      fwrite(fid,sprintf('%s :\t%g\t%g\t%g\n',abridgedTargetsName{i},abridgedTargetsNE(i,1),abridgedTargetsNE(i,2),abridgedTargetsNE(i,3)));
      labels = sprintf('%s\t%s',labels,abridgedTargetsName{i});
   end
   fwrite(fid,sprintf('Labels\n'));
   fwrite(fid,sprintf('%s\n',labels(2:end)));
   fclose('all');
   %% assemble and write eeglab chanlocs
   easyCapFCmmRigid = porcupineSequentialRotateAboutAxis(easyCapFCmm,-1,FIDSmm(6,:),alpha,beta,gamma);
   easyCapFCmmRigidO = easyCapFCmmRigid(1:64,:) - repmat(FIDSmm(6,:),[64,1]);
   fclose('all');
   fid = fopen(strcat(ntffp,'CHANLOCS_BCTMS64X13M43V3.xyz'),'w');
   for i=1:length(easyCapFCmmRigidO)
      fwrite(fid,sprintf('%02u %0.5f %0.5f %0.5f %s\n',i,-easyCapFCmmRigidO(i,2),easyCapFCmmRigidO(i,1),easyCapFCmmRigidO(i,3),easyCapLabels{i}));
   end
   fclose('all');
   % assemble and write brainstorm fiducials (milimeters), chanpos (centimeters), and scouts
   % convert to 256, 1mm neuroLRAPIS isospace
   ACvx = vtm(vtm(:,1) == 130,2:end);
   PCvx = vtm(vtm(:,1) == 140,2:end);
   IHvx = vtm(vtm(:,1) == 150,2:end);
   brainstormFiducialsNames = {'NAS','LPA','RPA','AC','PC','IH'};
   brainstormFiducialsVX = [NASvx;LPAvx;RPAvx;ACvx;PCvx;IHvx];
   brainstormFiducialsMM = porcupineCoordXfmNeuroVoxToBrainstormMM([pixdim1,pixdim2,pixdim3],[size(t,1),size(t,2),size(t,3)],brainstormFiducialsVX);
   brainstormFiducialsCM = brainstormFiducialsMM ./ 10;
   brainstormEasyCapMM = porcupineCoordXfmNeuroVoxToBrainstormMM([pixdim1,pixdim2,pixdim3],[size(t,1),size(t,2),size(t,3)],easyCapFCvx(1:64,:));
   brainstormEasyCapCM = brainstormEasyCapMM ./ 10;
   fid = fopen(strcat(ntffp,'BRAINSTORM_FIDUCIALS.m'),'w');
   for i=1:size(brainstormFiducialsMM,1)
      fwrite(fid,sprintf('%s = [%g,%g,%g];\n',brainstormFiducialsNames{i},brainstormFiducialsMM(i,1),brainstormFiducialsMM(i,2),brainstormFiducialsMM(i,3)));
   end
   fclose('all');
   fid = fopen(strcat(ntffp,'BRAINSTORM_CHANPOS_BCTMS64X13M43V3.txt'),'w');
   for i=1:3
      fwrite(fid,sprintf('%s\t%g\t%g\t%g\n',brainstormFiducialsNames{i},brainstormFiducialsCM(i,1),brainstormFiducialsCM(i,2),brainstormFiducialsCM(i,3)));
   end
   for i=1:64
      fwrite(fid,sprintf('%s\t%g\t%g\t%g\n',easyCapLabels{i},brainstormEasyCapCM(i,1),brainstormEasyCapCM(i,2),brainstormEasyCapCM(i,3)));
   end
   fclose('all');
   % assemble and write brainstorm scout seeds
   brainstormScoutSeeds = sortrows(masterTargets,5);
   j=1;
   for i=1:size(brainstormScoutSeeds,1)
      if ~isnan(brainstormScoutSeeds{i,5})
         brainstormScoutSeedsLabels{j,1} = brainstormScoutSeeds{i,2};
         brainstormScoutSeedsVX(j,1:3) = vtm(vtm(:,1) == brainstormScoutSeeds{i,1},2:end);
         brainstormScoutSeedsID(j,1) = brainstormScoutSeeds{i,1};
         j=j+1;
      end
   end
   brainstormScoutSeedsMM = porcupineCoordXfmNeuroVoxToBrainstormMM([pixdim1,pixdim2,pixdim3],[size(t,1),size(t,2),size(t,3)],brainstormScoutSeedsVX);
   save(strcat(ntffp,'BSS.mat'),'brainstormScoutSeedsLabels','brainstormScoutSeedsID','brainstormScoutSeedsMM')
   %% write images
   markerValInc = 10^floor(log10(max(t(:))));
   markerValBase = ceil(max(t(:)) / markerValInc) * markerValInc;
   tCapOverlay = t;
   tCap = zeros(size(t));
   for i=1:64
      x = easyCapFCvx(i,1);
      y = easyCapFCvx(i,2);
      z = easyCapFCvx(i,3);
      xx = max(1,x-1):min(size(t,1),x+1);
      yy = max(1,y-1):min(size(t,2),y+1);
      zz = max(1,z-1):min(size(t,3),z+1);
      tCapOverlay(xx,yy,zz) = markerValBase + i;
      tCap(x,y,z) = i;
   end
   tScoutsOverlay = t;
   tScouts = zeros(size(t));
   for i=1:size(brainstormScoutSeedsVX,1)
      x = brainstormScoutSeedsVX(i,1);
      y = brainstormScoutSeedsVX(i,2);
      z = brainstormScoutSeedsVX(i,3);
      xx = max(1,x-1):min(size(t,1),x+1);
      yy = max(1,y-1):min(size(t,2),y+1);
      zz = max(1,z-1):min(size(t,3),z+1);
      tScoutsOverlay(xx,yy,zz) = markerValBase + brainstormScoutSeedsID(i);
      tScouts(x,y,z) = brainstormScoutSeedsID(i);
   end
   porcupineVoxToNiiStan(tCapOverlay,sdn,strcat(ntffp,'VECO_BCTMS64X13M43V3'),fsldir);
   porcupineVoxToNiiStan(tCap,sdn,strcat(ntffp,'VEC_BCTMS64X13M43V3'),fsldir);
   porcupineVoxToNiiStan(tScoutsOverlay,sdn,strcat(ntffp,'BSSO'),fsldir);
   porcupineVoxToNiiStan(tScouts,sdn,strcat(ntffp,'BSS'),fsldir);
   tCapOverlay = t;
   for i=1:length(meridianCapFCvx)
      x = meridianCapFCvx(i,1);
      y = meridianCapFCvx(i,2);
      z = meridianCapFCvx(i,3);
      xx = max(1,x-1):min(size(t,1),x+1);
      yy = max(1,y-1):min(size(t,2),y+1);
      zz = max(1,z-1):min(size(t,3),z+1);
      tCapOverlay(xx,yy,zz) = markerValBase + NAStoINIcm;
   end
   for i=1:length(equatorCapFCvx)
      x = equatorCapFCvx(i,1);
      y = equatorCapFCvx(i,2);
      z = equatorCapFCvx(i,3);
      xx = max(1,x-1):min(size(t,1),x+1);
      yy = max(1,y-1):min(size(t,2),y+1);
      zz = max(1,z-1):min(size(t,3),z+1);
      tCapOverlay(xx,yy,zz) = markerValBase + LPAtoRPAcm;
   end
   for i=1:length(circumferenceCapFCvx)
      x = circumferenceCapFCvx(i,1);
      y = circumferenceCapFCvx(i,2);
      z = circumferenceCapFCvx(i,3);
      xx = max(1,x-1):min(size(t,1),x+1);
      yy = max(1,y-1):min(size(t,2),y+1);
      zz = max(1,z-1):min(size(t,3),z+1);
      tCapOverlay(xx,yy,zz) = markerValBase + CIRCcm;
   end
   porcupineVoxToNiiStan(tCapOverlay,sdn,strcat(ntffp,'VHM'),fsldir);
   status = true;
end