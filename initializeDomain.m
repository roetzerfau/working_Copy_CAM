 function [g, bulkVector,bulkTypeVector,chargeInfo, concAgent, edgeChargeVector, particleTypeVector, POMVector, ...
     POMconcVector, concPOMAgent, POMagentAge, solidParticleList, POMParticleList, randomPOMparticles, ...
     randomPOMparticlesSizes, POMsolidEdgeList] = initializeDomain(geoInputTxt, particlesInputTxt, randomPOMinputTxt, startConcPOM)  

% read geometry and particle indices from txt files
fileID = fopen(geoInputTxt,'r');
formatSpec = '%f';
geo = fscanf(fileID,formatSpec);
fclose(fileID);

fileID = fopen(particlesInputTxt,'r');
particleTypes = fscanf(fileID,formatSpec);
fclose(fileID);

numSquares = size(geo,1);
sizeStruc = sqrt(size(geo,1));

% create domain
intBound    = sizeStruc * ones(sizeStruc+1,1);
upBound     = intBound;

g = createDomainFolded(sizeStruc, sizeStruc, sizeStruc, 0, intBound, upBound);

% initialize vectors
bulkVector      = zeros(numSquares, 1);
bulkTypeVector  = zeros(numSquares, 1);
concAgent       = 0*ones( g.numE , 1 );   
edgeChargeVector = 0*ones( g.numCE , 1);    
chargeInfo      = zeros(numSquares, 2);
concPOMAgent = 0*ones( g.numCE , 1 ); 
POMagentAge = 0*ones( g.numCE , 1 ); 
POMVector = zeros(numSquares,1);


% reshape vectos, since different order of indices are used in input
sizeStruc = sqrt(size(geo,1));
geo = reshape(geo,[sizeStruc,sizeStruc]);
geo = rot90(geo,3);
geo = reshape(geo,[sizeStruc*sizeStruc,1]);

particleTypes = reshape(particleTypes,[sizeStruc,sizeStruc]);
particleTypes = rot90(particleTypes,3);
particleTypes = reshape(particleTypes,[sizeStruc*sizeStruc,1]);

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
POMconcVector = startConcPOM * POMVector;

% create list of unseparable solid particles directly from particleTypeVector 
solidParticleList = createSolidParticleList(particleTypeVector);

% create List of POM particles by finding connected components in POMVector
[POMParticleList, ~] = createPOMParticleList(POMVector);

% read POM shapes from txt file
POMshapesVector = readTxt(randomPOMinputTxt);

% create list of POM particles for new input by finding connected
% components in POMshapesVector. randomPOMparticlesSizes refers to the
% maximum size of the bounding box
[randomPOMparticles, randomPOMparticlesSizes] = createPOMParticleList(POMshapesVector);

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
[bulkVector, POMVector, POMParticleList, POMsolidEdgeList] = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);


end
