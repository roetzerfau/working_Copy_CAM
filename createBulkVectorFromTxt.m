 function [bulkVector,bulkTypeVector,chargeInfo, concAgent, edgeChargeVector, particleTypeVector, POMVector] = ...
    createBulkVectorFromTxt( g  , concAgent , NZd , edgeChargeVector )  

% rng('shuffle');
% rng(2);
numSquares      = (g.NX)^2; 
bulkVector      = zeros(numSquares, 1);
bulkTypeVector  = zeros(numSquares, 1);
particleTypeVector  = zeros(numSquares, 1);
chargeInfo      = zeros(numSquares, 2);

% fileID = fopen('Input/geo256ExKlebepunkteNew_surf15_new.txt','r');
fileID = fopen('Input/TUMgeoNew.txt','r');
formatSpec = '%f';
geo = fscanf(fileID,formatSpec);

% fileID = fopen('Input/geo256ExKlebepunkteNew_surf15_new_particles.txt','r');
fileID = fopen('Input/TUMparticlesNew.txt','r');
particleTypes = fscanf(fileID,formatSpec);

sizeStruc = sqrt(size(geo,1));
geo = reshape(geo,[sizeStruc,sizeStruc]);
geo = rot90(geo,3);
geo = reshape(geo,[sizeStruc*sizeStruc,1]);

particleTypes = reshape(particleTypes,[sizeStruc,sizeStruc]);
particleTypes = rot90(particleTypes,3);
particleTypes = reshape(particleTypes,[sizeStruc*sizeStruc,1]);

% bulkVector(strucMatrix == 1) = 1;
% bulkTypeVector(strucMatrix == 1) = -1;

bulkVector( geo == 1 ) = 1;
bulkVector( geo == 2 ) = 2;
% zu ändern! geo == 99 sollte nicht solid sein
bulkVector( geo == 99 ) = 99;

edgeChargeVector = setEdgeChargeVector(g, bulkVector, edgeChargeVector);

bulkVector( geo == 99 ) = 1;
bulkVector( geo == 2 ) = 1;

bulkTypeVector( geo == 1 ) = -1;
bulkTypeVector( geo == 99 ) = -1;
bulkTypeVector( geo == 2 ) = -1;

particleTypeVector = particleTypes;
POMVector = zeros(sizeStruc*sizeStruc,1);
POMVector( geo == 2 ) = 1;
% bulkVector = zeros(numSquares, 1);


% size(find(geo == 1))
% size(find(particleTypes > 0))
end
