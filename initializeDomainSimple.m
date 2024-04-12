 function [g, bulkVector,bulkTypeVector, concAgent, edgeChargeVector, reactiveSurfaceVector, particleTypeVector,...
     POMVector, POMconcVector, POMageVector, concPOMAgent, POMagentAge, solidParticleList, POMParticleList] ...
     = initializeDomainSimple(geo, porosity)  

sizeStruc = size(geo,1);
numSquares = sizeStruc*sizeStruc;

% create domain
intBound    = sizeStruc * ones(sizeStruc+1,1);
upBound     = intBound;

g = createDomainFolded(sizeStruc, sizeStruc, sizeStruc, 0, intBound, upBound);

% initialize vectorss
bulkVector      = zeros(numSquares, 1);
bulkTypeVector  = zeros(numSquares, 1);
concAgent       = 0*ones( g.numE , 1 );   
edgeChargeVector = 0*ones( g.numCE , 1);    
reactiveSurfaceVector = 0*ones( g.numCE , 1);  
concPOMAgent = 0*ones( g.numCE , 1 ); 
POMagentAge = 0*ones( g.numCE , 1 ); 
POMVector = zeros(numSquares,1);
POMageVector = zeros(numSquares,1);
particleTypes      = zeros(numSquares, 1);

numSolidPixels = round(porosity*numSquares);
solidInds = randperm(numSquares,numSolidPixels);

geo(solidInds) = 1;
for i = 1 : length(solidInds)
    particleTypes(solidInds(i)) = i;
end

% solidParticleList = cell(1,numSolidPixels);
%     
% for i = 1 : numSolidPixels
%     solidParticleList{i} = solidInds(i);
% end

% reshape vectos, since different order of indices are used in input
% sizeStruc = sqrt(size(geo,1));
% geo = reshape(geo,[sizeStruc,sizeStruc]);
% geo = rot90(geo,3);
% geo = reshape(geo,[sizeStruc*sizeStruc,1]);

% particleTypes = reshape(particleTypes,[sizeStruc,sizeStruc]);
% particleTypes = rot90(particleTypes,3);
% particleTypes = reshape(particleTypes,[sizeStruc*sizeStruc,1]);

% set bulkVector = 1 for solid, = 2 for POM 
bulkVector( geo == 1 ) = 1;
bulkVector( geo == 2 ) = 2;
% change maybe, geo == 99 should not be solid, or at least be careful!
bulkVector( geo == 99 ) = 99;

% set attractive edges on solid - fluid/bio boundary where geo == 99 in
% input
edgeChargeVector = setEdgeChargeVector(g, bulkVector, edgeChargeVector);

% set bulkVector = 1 for solid and POM, bulkTypeVector = -1 for solid and
% POM
bulkVector( geo == 99 ) = 1;
bulkVector( geo == 2 ) = 1;

bulkTypeVector( geo == 1 ) = -1;
bulkTypeVector( geo == 99 ) = -1;
bulkTypeVector( geo == 2 ) = -1;

% set particleTypeVector directly from input
particleTypeVector = particleTypes;

% set POMVector = 1 for POM, 0 elsewhere
POMVector( geo == 2 ) = 1;

% initialise POMconc uniformly with startConcPOM
POMconcVector = POMVector;

% create list of unseparable solid particles directly from particleTypeVector 
solidParticleList = createSolidParticleList(particleTypeVector);

% create List of POM particles by finding connected components in POMVector
[POMParticleList, ~] = createPOMParticleList(POMVector);

end
