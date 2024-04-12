 function [g, bulkVector,bulkTypeVector, concAgent, edgeChargeVector, reactiveSurfaceVector, particleTypeVector, ...
     POMVector, POMconcVector, POMageVector, concPOMAgent, POMagentAge, solidParticleList, POMParticleList, randomPOMparticles, ...
     randomPOMparticlesSizes] = initializeDomainFromMat(inputMat, randomPOMinputShapes)  

% read geometry and particle indices from mat file
load(inputMat, 'bulkVector', 'bulkTypeVector', 'edgeChargeVector', 'particleList', ...
    'particleTypeVector', 'POMconcVector', 'POMParticleList', 'POMVector', ...
      'reactiveSurfaceVector');
% load(inputMat);

numSquares = size(bulkVector,1);
sizeStruc = sqrt(size(bulkVector,1));

% create domain
intBound    = sizeStruc * ones(sizeStruc+1,1);
upBound     = intBound;

g = createDomainFolded(sizeStruc, sizeStruc, sizeStruc, 0, intBound, upBound);

% initialize vectors
POMageVector = POMVector;
concAgent       = 0*ones( g.numE , 1 );   
concPOMAgent = 0*ones( g.numCE , 1 ); 
POMagentAge = 0*ones( g.numCE , 1 ); 

% create list of unseparable solid particles directly from particleTypeVector 
solidParticleList = particleList;
% POMsolidEdgeList = particleSurfaceEdgeList;

% read POM shapes from txt file
% POMshapesVector = readTxt(randomPOMinputTxt);

% create list of POM particles for new input by finding connected
% components in POMshapesVector. randomPOMparticlesSizes refers to the
% maximum size of the bounding box
% [randomPOMparticles, randomPOMparticlesSizes] = createPOMParticleList(POMshapesVector);
load(randomPOMinputShapes);

% calculate number of cells in each particle
randomPOMparticlesAreas = zeros(1, length(randomPOMparticles));
for i = 1 : length(randomPOMparticles)
    randomPOMparticlesAreas(i) = length(randomPOMparticles{i});
end
% reorder list of POM particles for input according to their descending
% area
[~, indOrder] = sort(randomPOMparticlesAreas);
randomPOMparticlesSizes = randomPOMparticlesSizes(indOrder);
randomPOMparticles = randomPOMparticles(indOrder);

% calculate POMsolidEdgeList, only for output of initial state
% [bulkVector, POMVector, POMParticleList, POMsolidEdgeList] = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);


end
