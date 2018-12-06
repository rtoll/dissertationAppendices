function kernelReduced = kernelPCAreductionUsingEEG(kernelUnconstrained,eegData)
   sourceSpace = kernelUnconstrained * eegData;
   transformationMatrix = zeros(size(sourceSpace,1)/3, size(sourceSpace,1));
   for i = 1: size(sourceSpace,1)/3
      sourceSpaceTemp = sourceSpace((i-1)*3+1:i*3,:);
      [eigenvectors, eigenvalues] = eig(sourceSpaceTemp * sourceSpaceTemp');
      [~, sortIndices] = sort(diag(eigenvalues),'descend');
      transformationMatrix(i,(i-1)*3+1:i*3) = eigenvectors(:,sortIndices(1))';
   end
   kernelReduced = transformationMatrix * kernelUnconstrained;
end