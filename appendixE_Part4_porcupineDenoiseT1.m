function [status] = porcupineDenoiseT1(inputT1basename,basedir,fsldir)
   % This function cleans a T1 image, removing the empty space noise
   % in the image, leaving a smooth scalp/skin boundary. It also removes
   % wrapped noses so long as the nose doesn't come into contact with the back
   % of head. Two .mat files are necessary
   % to run this function, c6_256 and s4_256.mat, they are adjacency matrices
   % The raw T1 must be <= than 256x256x256 voxels and already be in standard orientation (RLPAIS)
   status = false;
   setenv('FSLDIR',fsldir);
   setenv('FSLOUTPUTTYPE', 'NIFTI');
   % define and create directories
   strawFolder = strcat(basedir,'/Processing/Stage2_StructuralRegistrationAndWarps/',inputT1basename,'_STRAW');
   system(sprintf('mkdir -p %s',strawFolder));
   % define filenames
   inputOVFFullPath = strcat(basedir,'/Processing/Stage1_StandardizeT1/',inputT1basename,'_OVF');
   outputStanDeNoisedFullPath = strcat(strawFolder,'/',inputT1basename,'_SDN'); % standardized and denoised
   % load the predefined adjacency matrices for 256^3 volumes
   load(strcat(basedir,'/Admin/c6_256.mat')) % volume
   load(strcat(basedir,'/Admin/s4_256.mat')) % slice
   % obtain dimensions
   [~,pixdim1] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',inputOVFFullPath,'pixdim1'));
   pixdim1 = str2double(pixdim1);
   [~,pixdim2] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',inputOVFFullPath,'pixdim2'));
   pixdim2 = str2double(pixdim2);
   [~,pixdim3] = system(sprintf('%s/%s %s %s',strcat(fsldir,'/bin'),'fslval',inputOVFFullPath,'pixdim3'));
   pixdim3 = str2double(pixdim3);
   % define thresholds
   threshBaseToPeakRatio = .1; % scalar, defines the height of the main histogram peak at which the peak is considered leveled off
   threshContrastRatio = 2; % scalar, minimum contrast
   threshGoodNeck = 25; % percentage, value of the sum of axial slices at which the neck is most prominent
   threshFaces = 1; % scalar, minimum number of adjacent positive faces needed to survive "sanding"
   threshClusterSizeMM2 = 400; % mm^2, planar clusters smaller than this are zeroed out
   threshMarginalVolumeChangeMM3 = 10; % mm^3, used to escape recursions when effective change is small
   threshClusterSizeX = round(threshClusterSizeMM2 / (pixdim2 * pixdim3)); % voxels conversion of threshClusterSizeMM2
   threshClusterSizeY = round(threshClusterSizeMM2 / (pixdim1 * pixdim3)); % voxels conversion of threshClusterSizeMM2
   threshClusterSizeZ = round(threshClusterSizeMM2 / (pixdim1 * pixdim2)); % voxels conversion of threshClusterSizeMM2
   threshMarginalVolumeChange = round(threshMarginalVolumeChangeMM3 / (pixdim1 * pixdim2 * pixdim3)); % voxel conversion of threshMarginalVolumeChangeMM3
   threshEscape  = .25; % ratio of cleaned image count of non zero voxels to original image, used to identify when sanding was too severe and revert to original
   % load original T1
   t = porcupineNiiToVoxStan(inputOVFFullPath,fsldir);
   % ensure all values are >= 0
   t = t - min(t(:));
   tb = false(size(t));
   % define standard value array
   stv = zeros(256,256,256);
   stv(1:size(t,1),1:size(t,2),1:size(t,3)) = t;
   % generate integer histogram by axial slice to standardize by peak value
   for z=1:size(t,3)
      tSlice = squeeze(t(:,:,z));
      if any(tSlice(:) ~= 0)
         [shix,shiy] = porcupineHistInt(tSlice(tSlice > 0));
         maxShiy = max(shiy);
         maxShix = shix(find(shiy == maxShiy,1,'first'));
         baseShiy = floor(maxShiy * threshBaseToPeakRatio);
         baseShix = shix(find(shiy <= baseShiy & shix > maxShix,1,'first'));
         if ~isempty(baseShix)
            tSlice(tSlice < baseShix) = baseShix;
            tSlice = tSlice ./ baseShix;
         else
            tSlice = ones(size(tSlice));
         end
      else
         tSlice = ones(size(tSlice));
      end
      tb(:,:,z) = tSlice >= threshContrastRatio;
   end
   % define standard binary array
   stb = false(256,256,256);
   stb(1:size(t,1),1:size(t,2),1:size(t,3)) = tb;
   % perfrom first contact indices sequence
   fc3 = porcupineFirstContactInclusive3D(stb);
   % seal all inner holes by only zeroing edge to fc3
   stb(:) = true;
   stb(fc3) = false;
   % find the first good axial slice and seal with square slab to first axial slice
   for z=1:size(t,3)
      tSlice = squeeze(tb(:,:,z));
      neckSlices(z,1) = sum(tSlice(:));
   end
   neckSlice = find(neckSlices >= prctile(neckSlices,threshGoodNeck),1,'first');
   neckSlice = min(neckSlice,size(t,3) * .25); % safety to prevent erroneous mid head cutoff as neck
   stb(:,:,1:neckSlice) = true;
   % iteratively remove voxels with fewer positive adjacent faces than thresh
   repeat = true;
   while repeat
      repeat = false;
      startSum = sum(stb(:));
      stb(sum(stb(c),2) <= threshFaces) = false;
      if startSum > sum(stb(:))
         repeat = true;
      end
   end
   % perfrom first contact indices sequence again 
   fc3 = porcupineFirstContactInclusive3D(stb);
   % seal all inner holes by only zeroing edge to fc3
   stb(:) = true;
   stb(fc3) = false;
   % remove neck slab and replicate first good neck slice + 1 to first axial slice
   for z=1:neckSlice
      stb(:,:,z) = stb(:,:,neckSlice + 1);
   end
   % zero boundary positive voxels to make tiny clusters easier to identify
   stb(any(~stb(c),2)) = false;
   % begin "sanding" of noisy artifact clusters smaller than thresh
   loops=1; % breakout in case of infinite sanding
   repeat1 = true;
   while repeat1 && loops <= 10
      repeat1 = false;
      startSum1 = sum(stb(:));
      % pixel cluster, pixels because this is done in 2d
      for stAxis=1:3
         for dim=1:256
            clc
            disp('PHASE 3 of 8: DENOISE T1')
            disp(sprintf('cluster iteration %u is %04.1f%% complete\n3 to 5 iterations are typical',loops,100*((stAxis-1)*256+dim)/(3*256)))
            if stAxis == 1
               stSlice = squeeze(stb(dim,:,:));
               threshClusterSize = threshClusterSizeX;
            elseif stAxis == 2
               stSlice = squeeze(stb(:,dim,:));
               threshClusterSize = threshClusterSizeY;
            else
               stSlice = squeeze(stb(:,:,dim));
               threshClusterSize = threshClusterSizeZ;
            end
            repeat2 = true;
            while repeat2
               repeat2 = false;
               startSum2 = sum(stSlice(:));
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % fc, 2d
               fc = []; % first contact
               for rowInd=1:256
                  straw = stSlice(rowInd,:);
                  hits = find(straw);
                  if ~isempty(hits);
                     if dim == 3
                        fc = [fc,sub2ind([256 256],rowInd,hits(1))];
                     end
                     fc = [fc,sub2ind([256 256],rowInd,hits(end))];
                  end
               end
               for colInd=1:256
                  straw = stSlice(:,colInd);
                  hits = find(straw);
                  if ~isempty(hits);
                     fc = [fc,sub2ind([256 256],hits(1),colInd)];
                     fc = [fc,sub2ind([256 256],hits(end),colInd)];
                  end
               end
               fc = unique(fc);
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               % cluster
               while ~isempty(fc)
                  tempCluster = fc(1);
                  s4cpi = stSlice(s4) .* s4; % connected positive indices
                  repeat3 = true;
                  while repeat3
                     repeat3 = false;
                     startLength3 = length(tempCluster);
                     cpi = s4cpi(tempCluster,:);
                     cpi = [tempCluster,cpi];
                     if size(cpi,2) > size(cpi,1)
                        cpi = cpi';
                     end
                     cpi = unique(cpi(cpi > 0));
                     tempCluster = cpi;
                     if startLength3 ~= length(tempCluster)
                        repeat3 = true;
                     end
                     if length(tempCluster) >= threshClusterSize
                        repeat3 = false;
                     end
                  end      
                  if length(tempCluster) < threshClusterSize
                     stSlice(tempCluster) = false;
                  end
                  fc = setdiff(fc,tempCluster);
               end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if startSum2 > sum(stSlice(:))
                  repeat2 = true;
               end
            end
            if stAxis == 1
               stb(dim,:,:) = stSlice;
            elseif stAxis == 2
               stb(:,dim,:) = stSlice;
            else
               stb(:,:,dim) = stSlice;
            end
         end
      end
      if startSum1 - sum(stb(:)) > threshMarginalVolumeChange
         repeat1 = true;
      end
      loops=loops+1;
   end
   % anneal gaps caused by excessive sanding
   for z=1:256
      stSlice = squeeze(stb(:,:,z));
      repeat = true;
      while repeat
         repeat = false;
         startSum = sum(stSlice(:));
         for x=1:256
            if any(stSlice(x,:))
               hits = find(stSlice(x,:));
               stSlice(x,hits(1):hits(end)) = true;
            end
         end
         for y=1:256
            if any(stSlice(:,y))
               hits = find(stSlice(:,y));
               stSlice(hits(1):hits(end),y) = true;
            end
         end
         if startSum < sum(stSlice(:));
            repeat = true;
         end
      end
      stb(:,:,z) = stSlice;
   end
   % take the 6 dof mean twice to smooth out jagged edges
   sts = single(stb);
   sts(:) = mean(sts(c),2);
   sts(:) = mean(sts(c),2); % intentional repeat
   % multiply original values by mask
   sts(sts > 0) = 1;
   stm = sts .* (stv + 1); % 1 is added to stv to prevent zero voxels within brain volume
   % perform fc3 sequence again 
   stb = stm > 0;
   fc3 = porcupineFirstContactInclusive3D(stb);
   stb = true(size(stm));
   stb(fc3) = false;
   sts = single(stb);
   sts = sts .* (stv + 1);
   % escape if sanding was too severe
   if sum(sts(:) > 0) < threshEscape * sum(stv(:) > 0)
      stm = stv;
   end
   % write out the final cleaned T1 as mat and image
   t = stm(1:size(t,1),1:size(t,2),1:size(t,3));
   save(strcat(outputStanDeNoisedFullPath,'.mat'),'t')
   porcupineVoxToNiiStan(t,inputOVFFullPath,outputStanDeNoisedFullPath,fsldir);
   status = true;
end