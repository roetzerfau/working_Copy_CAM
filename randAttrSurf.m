edgeChargeVector = 0*ones( g.numCE , 1);   

particleSurfaceHelper = zeros(length(particleSurfaceEdgeList));
for i = 1 : length(particleSurfaceEdgeList)
    particleSurfaceHelper(i) = length(particleSurfaceEdgeList{i});
end

numAttrEdges = 106;
for i = 1 : numAttrEdges
   particleCandidatesIndices = find(particleSurfaceHelper > 80);
   success = 0;
   while ~success
      particleInd = randi(length(particleCandidatesIndices));
      particleSurface = particleSurfaceEdgeList{particleCandidatesIndices(particleInd)};
      ind = randi(length(particleSurface));
      edgeCandidates = ind : (ind + 14);
      edgeCandidates(edgeCandidates> length(particleSurface)) = ...
          edgeCandidates(edgeCandidates> length(particleSurface)) - length(particleSurface);
      if sum(edgeChargeVector(particleSurfaceEdgeList{particleCandidatesIndices(particleInd)}(edgeCandidates))) == 0
          edgeChargeVector(particleSurfaceEdgeList{particleCandidatesIndices(particleInd)}(edgeCandidates)) = 1;
          success = 1;
      end
   end
end

for i = 1 : numAttrEdges
   particleCandidatesIndices = find(particleSurfaceHelper > 30);
   success = 0;
   while ~success
      particleInd = randi(length(particleCandidatesIndices));
      particleSurface = particleSurfaceEdgeList{particleCandidatesIndices(particleInd)};
      ind = randi(length(particleSurface));
      edgeCandidates = ind : (ind + 9);
      edgeCandidates(edgeCandidates> length(particleSurface)) = ...
          edgeCandidates(edgeCandidates> length(particleSurface)) - length(particleSurface);
      if sum(edgeChargeVector(particleSurfaceEdgeList{particleCandidatesIndices(particleInd)}(edgeCandidates))) == 0
          edgeChargeVector(particleSurfaceEdgeList{particleCandidatesIndices(particleInd)}(edgeCandidates)) = 1;
          success = 1;
      end
   end
end

% for i = 1 : numAttrEdges
%    particleCandidatesIndices = find(particleSurfaceHelper > 20);
%    success = 0;
%    while ~success
%       particleInd = randi(length(particleCandidatesIndices));
%       particleSurface = particleSurfaceEdgeList{particleCandidatesIndices(particleInd)};
%       ind = randi(length(particleSurface));
%       edgeCandidates = ind : (ind + 4);
%       edgeCandidates(edgeCandidates> length(particleSurface)) = ...
%           edgeCandidates(edgeCandidates> length(particleSurface)) - length(particleSurface);
%       if sum(edgeChargeVector(particleSurfaceEdgeList{particleCandidatesIndices(particleInd)}(edgeCandidates))) == 0
%           edgeChargeVector(particleSurfaceEdgeList{particleCandidatesIndices(particleInd)}(edgeCandidates)) = 1;
%           success = 1;
%       end
%    end
% end


bulkVector = particleTypeVector;
bulkVector(bulkVector>0)=1;
visualizeDataSub(g, bulkVector, 'bulk', 'test', 0);
visualizeDataEdges(g, edgeChargeVector, 'conc', 'test2', 0, 2);