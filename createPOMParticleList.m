function [particleList, particleSizesHelper] = createPOMParticleList (POMVector) 

if size(POMVector,1) ~= size(POMVector,2)
    sizeStruc = sqrt(length(POMVector));
    POMVector = reshape(POMVector,[sizeStruc,sizeStruc]);
end

% POMVector = rot90(POMVector,1);

cc = bwconncomp(POMVector,4);

POMlabels = bwlabel(POMVector, 4);
POMboundingBoxes = regionprops(POMlabels,'BoundingBox');
POMboundingBoxes = struct2cell(POMboundingBoxes);

particleSizesHelper = zeros(length(POMboundingBoxes),1);
for i = 1 : length(POMboundingBoxes)
    particleSizesHelper(i) = max( POMboundingBoxes{i}(3),  POMboundingBoxes{i}(4) );
end

particleList = cc.PixelIdxList;

for i = 1 : length(particleList)
   particleList{i} = particleList{i}'; 
end
%     numSolidParticles = max(max(particleTypeVector));
%     particleList = cell(1,numSolidParticles);
%     
%     for i = 1 : numSolidParticlesPOMconc_temp
%         indSolidParticle = find(particleTypeVector == i);
%         particleList{i} = indSolidParticle';
%     end
end
