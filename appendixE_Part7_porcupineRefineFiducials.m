function [status] = porcupineRefineFiducials(inputT1basename,basedir,fsldir)
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
   % load planar pixel adjacency matrix
   load(strcat(basedir,'/Admin/s4_256.mat')) % slice
   % obtain T1 dimensions
   [~,pixdim1] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',sdn,'pixdim1'));
   pixdim1 = str2double(pixdim1);
   [~,pixdim2] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',sdn,'pixdim2'));
   pixdim2 = str2double(pixdim2);
   [~,pixdim3] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',sdn,'pixdim3'));
   pixdim3 = str2double(pixdim3);
   % load sdn
   t = porcupineNiiToVoxStan(sdn,fsldir);
   % load sdn into 256 space
   st = zeros(256,256,256);
   st(1:size(t,1),1:size(t,2),1:size(t,3)) = t;
   stVal = st;
   stBin = logical(st);
   % calculate sanding thresholds
   stt = st(:);
   stt = stt(stt > 0);
   threshNAS = prctile(stt,50);
   threshLPAandRPA = prctile(stt,25);
   % load voxel targets matrix
   load(strcat(ntffp,'VTM.mat'))
   % load initial fiducial estimates (lower case are estimates, upper case are refined)
   nas = vtm(vtm(:,1) == 110,2:end);
   ini = vtm(vtm(:,1) == 160,2:end);
   lpa = vtm(vtm(:,1) == 121,2:end);
   rpa = vtm(vtm(:,1) == 122,2:end);
   % define distance error threshold
   threshDist = 10; % mm, if the refined position is farther than threshDist from estimated position,
                    % assume error occurred and use estimate as the refined value
   % load fixed PC location
   PC = vtm(vtm(:,1) == 140,2:end);
   % nasion refinement
   threshSmallClusterMM2 = 9; % mm^2, planar cluster area
   threshSmallCluster = round(threshSmallClusterMM2 / (pixdim1 * pixdim3));
   nudgeYaMM = 1; % mm, once nasion gap found, nudge anterior to skin surface
   sclXMM = 20; % mm, search cube length, x dim
   sclYMM = 20; % mm, search cube length, y dim
   sclZMM = 20; % mm, search cube length, z dim
   nudgeYa = round(nudgeYaMM / pixdim2);
   sclX = round(sclXMM / pixdim1);
   sclY = round(sclYMM / pixdim2);
   sclZ = round(sclZMM / pixdim3);
   st = zeros(size(st));
   nasXmin = max(1,nas(1) - sclX);
   nasXmax = min(256,nas(1) + sclX);
   nasYmin = max(1,nas(2) - sclY);
   nasYmax = min(256,nas(2) + sclY);
   nasZmin = max(1,nas(3) - sclZ);
   nasZmax = min(256,nas(3) + sclZ);
   nasXrange = nasXmin:nasXmax;
   nasYrange = nasYmin:nasYmax;
   nasZrange = nasZmin:nasZmax;
   st(nasXrange,nasYrange,nasZrange) = stVal(nasXrange,nasYrange,nasZrange);
   st(st < threshNAS) = 0;
   st = logical(st);
   for y=nasYrange
      stSlice = squeeze(st(:,y,:));
      % cluster
      nzi = find(stSlice);
      clusters = zeros(length(nzi),1);
      while ~isempty(nzi)
         tc = nzi(1);
         repeat = true;
         while repeat
            repeat = false;
            startLength = length(tc);
            cpi = s4(tc,:) .* stSlice(s4(tc,:));
            cpi = cpi(cpi > 0);
            tc = unique([tc(:);cpi(:)]);
            if startLength < length(tc)
               repeat = true;
            end      
         end
         clusters(1:length(tc),end) = tc;
         clusters(:,end+1) = 0;
         nzi = setdiff(nzi,tc);
      end
      clusters = clusters(:,1:end-1);
      keepIndices = [];
      for i=1:size(clusters,2)
         tc = clusters(:,i);
         tc = tc(tc > 0);
         if length(tc) >= threshSmallCluster
            keepIndices = [keepIndices;tc];
         end
      end
      stSlice = false(256,256);
      stSlice(keepIndices) = true;
      st(:,y,:) = stSlice;
   end
   i=1;
   for x=nasXrange
      for z=nasZrange
         y = find(st(x,:,z),1,'last');
         if ~isempty(y)
            st(x,nasYmin:y,z) = true;
            fc(i,:) = [x,y,z];
            i=i+1;
         end
      end
   end
   firstGapY = [];
   advance = true;
   y=nasYmin;
   while advance && y <= nasYmax
      stSlice = squeeze(st(:,y,:));
      for i=1:length(nasZrange)
         zeroZs(i,1) = sum(all(~stSlice(:,nasZrange(i))));
      end
      demark = round((2*sclX+1)/3);
      if sum(zeroZs(1:demark) == 0) > demark/2 && sum(zeroZs(demark+1:2*demark) == 1) > demark/2 && sum(zeroZs(2*demark+1:end) == 0) > demark/2
         firstGapY = y;
         advance = false;
      end
      y=y+1;
   end
   if ~isempty(firstGapY)
      y = firstGapY;
   else
      y = nas(2);
      sprintf('WARNING: did not locate a nasion gap')
   end
   z = round(nasZmin + median(find(zeroZs)));
   if isnan(z)
      z = nas(3);
   end
   x = round(median(find(squeeze(st(:,y-1,z)))));
   if isnan(x)
      x = nas(1);
   end
   st = stBin;
   if st(x,y+nudgeYa,z)
      NAS = [x,y+nudgeYa,z];
   elseif st(x,y+round(nudgeYa/2),z)
      NAS = [x,y+round(nudgeYa/2),z];
   elseif st(x,y,z)      
      NAS = [x,y,z];
   else
      NAS = [x,y-1,z];
   end
   distNASnas = sqrt((nas(1) * pixdim1 - NAS(1) * pixdim1)^2 + (nas(2) * pixdim2 - NAS(2) * pixdim2)^2 + (nas(3) * pixdim3 - NAS(3) * pixdim3)^2);
   if distNASnas > threshDist
      NAS = nas;
      sprintf('WARNING: refined - estimate distance for nasion exceeds threshold, using estimate value')
   end
   % append NAS to vtm
   vtm(end+1,:) = [210,NAS];
   % inion refinement
   nudgeZiMM = 2; % mm, once inion found, nudge inferior to below cavity's center of mass
   nudgeZi = round(nudgeZiMM / pixdim3);
   st = stBin;
   m = (NAS(2) - PC(2)) / (NAS(1) - PC(1) + .001); % .001 added to prevent infinite slope
   b = NAS(2) - (m * NAS(1));
   y = 1;
   x = round((y - b) / m);
   z = ini(3) - nudgeZi;
   if z < 1
      z = 1;
   end
   if st(x,1,z)
      % if the first y slice is a hit, then assume bad cleaning and use
      % stVal at threshNAS
      st = stVal;
      st(st < threshNAS) = 0;
      st = logical(st);
   end
   while ~st(x,y,z)
      y=y+1;
      x = round((y - b) / m);
   end
   INI = [x,y,z];
   % INI is excluded from threshDist because it relies on NAS and is scooted back from within skull
   % to skin surface
   % append INI to vtm
   vtm(end+1,:) = [260,INI];
   % left and right preauricular refinement
   srlXMM = 20; % mm, search rectangle length, x dim
   srlYaMM = 10; % mm, search rectangle length, y dim, anterior
   srlYpMM = 20; % mm, search rectangle length, y dim, posterior
   srlZMM = 20; % mm, search rectangle length, z dim
   nudgeYaMM = 2; % mm, once tragus found, nudge anterior to center mass
   smoothingLengthMM = 3; % mm, distance to smooth ear canal sum
   srlX = round(srlXMM / pixdim1);
   srlYa = round(srlYaMM / pixdim2);
   srlYp = round(srlYpMM / pixdim2);
   srlZ = round(srlXMM / pixdim3);
   if min(lpa(3),rpa(3)) - srlZ < 1
      srlZ = min(lpa(3),rpa(3)) - 1;
   end
   nudgeYa = round(nudgeYaMM / pixdim2);
   smoothingLength = round(smoothingLengthMM / pixdim3);
   indexx = 1:256;
   st = stVal;
   st(st < threshLPAandRPA) = 0;
   st = logical(st);
   for ear=1:2
      rectMask = false(256,256);
      bi = [];
      if ear == 1
         % left
         y = lpa(2);
         z0 = lpa(3);
         x = find(st(:,y,z0),1,'first');
         y = find(~st(x,:,z0) & indexx < y,1,'last') + 1;      
         rectMask(x:x + srlX,y - srlYp:y + srlYa) = true;
      else
         % right
         y = rpa(2);
         z0 = rpa(3);
         x = find(st(:,y,z0),1,'last');
         y = find(~st(x,:,z0) & indexx < y,1,'last') + 1;
         rectMask(x - srlX:end,y - srlYp:y + srlYa) = true;
      end
      for y=1:256
         bi = [bi,sub2ind([256 256],x,y)];
      end
      axialEarCanalClusterSum = zeros(256,1);
      for zz=-srlZ:srlZ
         z = z0 + zz;
         stSlice = squeeze(~st(:,:,z));
         stSlice = stSlice .* rectMask;
         stSlice(~any(stSlice(s4),2)) = false;
         nzi = find(stSlice);
         clusters = zeros(length(nzi),1);
         while ~isempty(nzi)
            tc = nzi(1);
            repeat = true;
            while repeat
               repeat = false;
               startLength = length(tc);
               cpi = s4(tc,:) .* stSlice(s4(tc,:));
               cpi = cpi(cpi > 0);
               tc = unique([tc(:);cpi(:)]);
               if startLength < length(tc)
                  repeat = true;
               end      
            end
            clusters(1:length(tc),end) = tc;
            clusters(:,end+1) = 0;
            nzi = setdiff(nzi,tc);
         end
         clusters = clusters(:,1:end-1);
         clustersBi = zeros(size(clusters,1),1);
         for i=1:size(clusters,2)
            if ~isempty(intersect(clusters(:,i),bi))
               clustersBi(:,end) = clusters(:,i);
               clustersBi(:,end+1) = 0;
            end
         end
         clustersBi = clustersBi(:,1:end-1);
         if ~isempty(clustersBi)
            clustersBiSum = sum(clustersBi > 0);
            clustersBi = clustersBi(:,find(clustersBiSum == max(clustersBiSum),1,'first'));
            clustersBi = clustersBi(clustersBi > 0);
            axialEarCanalClusterSum(z,1) = length(clustersBi);
         end
      end
      axialEarCanalClusterSumSmoothed = zeros(size(axialEarCanalClusterSum));
      for i=1:256
         if i > smoothingLength && i <= 256 - smoothingLength
            axialEarCanalClusterSumSmoothed(i,1) = mean(axialEarCanalClusterSum(i-smoothingLength:i+smoothingLength));
         end
      end
      if ear == 1
         y = lpa(2);
         z = find(axialEarCanalClusterSumSmoothed == max(axialEarCanalClusterSumSmoothed),1,'first');
         x = find(st(:,y,z),1,'first');
         y = find(~st(x,:,z) & indexx < y,1,'last') + 1;
         if st(x,y + nudgeYa,z)
            y = y + nudgeYa;
            x = find(st(:,y,z),1,'first');
         elseif st(x,round(y+nudgeYa/2),z)
            y = round(y + nudgeYa / 2);
            x = find(st(:,y,z),1,'first');
         else
            x = find(st(:,y,z),1,'first');
         end
         LPA = [x y z];
         distLPAlpa = sqrt((lpa(1) * pixdim1 - LPA(1) * pixdim1)^2 + (lpa(2) * pixdim2 - LPA(2) * pixdim2)^2 + (lpa(3) * pixdim3 - LPA(3) * pixdim3)^2);
         if distLPAlpa > threshDist
            LPA = lpa;
            sprintf('WARNING: refined - estimate distance for LPA exceeds threshold, using estimate value')
         end
      else
         y = rpa(2);
         z = find(axialEarCanalClusterSumSmoothed == max(axialEarCanalClusterSumSmoothed),1,'first');
         x = find(st(:,y,z),1,'last');
         y = find(~st(x,:,z) & indexx < y,1,'last') + 1;
         if st(x,y + nudgeYa,z)
            y = y + nudgeYa;
            x = find(st(:,y,z),1,'last');
         elseif st(x,round(y + nudgeYa / 2),z)
            y = round(y + nudgeYa / 2);
            x = find(st(:,y,z),1,'last');
         else
            x = find(st(:,y,z),1,'last');
         end
         RPA = [x y z];
         distRPArpa = sqrt((rpa(1) * pixdim1 - RPA(1) * pixdim1)^2 + (rpa(2) * pixdim2 - RPA(2) * pixdim2)^2 + (rpa(3) * pixdim3 - RPA(3) * pixdim3)^2);
         if distRPArpa > threshDist
            RPA = rpa;
            sprintf('WARNING: refined - estimate distance for RPA exceeds threshold, using estimate value')
         end
      end
   end
   % append LPA and RPA to vtm
   vtm(end+1,:) = [221,LPA];
   vtm(end+1,:) = [222,RPA];
   % re-sort vtm
   vtm = sortrows(vtm,1);
   vtmRad = [size(t,1) - vtm(:,1) + 1,vtm(:,2),vtm(:,3)];
   % overwrite previous vtms with updated versions
   save(strcat(ntffp,'VTM'),'vtm') % voxel targets matrix
   save(strcat(ntffp,'VTMrad'),'vtmRad') % voxel targets matrix radiological
   status = true;
end