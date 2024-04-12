 function [g, bulkVector, bulkTypeVector, reactiveSurfaceVector, particleTypeVector,...
    solidParticleList] = initializeDomainRandom(N, porosity)

% read geometry and particle indices from mat file

sizeStruc = N;
numSquares = N*N;

% create domain
intBound    = sizeStruc * ones(sizeStruc+1,1);
upBound     = intBound;

g = createDomainFolded(sizeStruc, sizeStruc, sizeStruc, 0, intBound, upBound);

% initialize vectors
bulkVector      = zeros(numSquares, 1);
bulkTypeVector  = zeros(numSquares, 1);
particleTypeVector  = zeros(numSquares, 1);
reactiveSurfaceVector = 0*ones( g.numCE , 1);    

numSolidPixels = round(porosity*numSquares);
solidInds = randperm(numSquares,numSolidPixels);

bulkVector(solidInds) = 1;
bulkTypeVector(solidInds) = -1;
particleTypeVector(solidInds) = 1;

solidParticleList = cell(1,numSolidPixels);
    
for i = 1 : numSolidPixels
    solidParticleList{i} = solidInds(i);
end

end
