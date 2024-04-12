 function [bulkVector,bulkTypeVector,chargeInfo,particle1List, particle2List, particle3List, particle4List,...
    particle5List, particle6List, concAgent, edgeChargeVector, particleTypeVector] = ...
    createBulkVector( g , porosity , concAgent , NZd , amountOfParticles , chargeAtEdges, edgeChargeVector )  

rng('shuffle');
% rng(2);
numSquares      = (g.NX)^2; 
bulkVector      = zeros(numSquares, 1);
bulkTypeVector  = zeros(numSquares, 1);
particleTypeVector  = zeros(numSquares, 1);
chargeInfo      = zeros(numSquares, 2);
numFreeSquares    = floor(numSquares*(1-porosity)); 

% particle1 1x17  (Goethit)
% particle2 2x34  (Goethit)
% particle3 2x6   (Illite)
% particle4 2x30  (Illite)
% particle5 8x100 (Illite)
% particle6 1x1

numParticle1 = floor( (numFreeSquares*amountOfParticles(1)/17) );
numParticle2 = floor( (numFreeSquares*amountOfParticles(2)/68) );
numParticle3 = floor( (numFreeSquares*amountOfParticles(3)/12) );
numParticle4 = floor( (numFreeSquares*amountOfParticles(4)/60) );
numParticle5 = floor( (numFreeSquares*amountOfParticles(5)/800));
numParticle6 = floor( (numFreeSquares*amountOfParticles(6)));

particle1List = zeros(numParticle1,17);
particle2List = zeros(numParticle2,68);
particle3List = zeros(numParticle3,12);
particle4List = zeros(numParticle4,60);
particle5List = zeros(numParticle5,800);
particle6List = zeros(numParticle6,1);

surfaceParticle1 = 36;
surfaceParticle2 = 72;
surfaceParticle3 = 16;
surfaceParticle4 = 64;
surfaceParticle5 = 216;
surfaceParticle6 = 4;

stencilSize = [10, 18, 4, 16, 53, 3];

for i = [5,2,4,1,3,6] % order, in which the particles are placed (largest to smallest particle)

if i == 1
%% place particle 1
fprintf('Start placing particle 1 (1x17): \n')
sumParticle = 0;
numVertical = floor(numParticle1/2);
numHorizontal = numParticle1-numVertical;
verticalInds = findParticleInds(17,1); % local index in stencil
horizontalInds = findParticleInds(1,17);

%place vertical particles
while sumParticle < numVertical
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(verticalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = 1;
        particleTypeVector(particleInds) = 1;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(1,:),17,1);  %[kurze, lange] 
        particle1List(sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 1 1], 17, 1);     
    end
end

sumParticle = 0;    
%place horizontal particles
while sumParticle < numHorizontal
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(horizontalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = 1;
        particleTypeVector(particleInds) = 1;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(1,:),17,1);
        particle1List(numVertical+sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 1 1], 17, 1);
    end
end    
fprintf('%d particle 1 placed \n',numVertical+numHorizontal)   
end
if i == 2
%% place particle 2
fprintf('Start placing particle 2 (2x34): \n')
sumParticle = 0;
numVertical = floor(numParticle2/2);
numHorizontal = numParticle2-numVertical;
verticalInds = findParticleInds(34,2); % local index in stencil
horizontalInds = findParticleInds(2,34);

%place vertical particles
while sumParticle < numVertical
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(verticalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = 1;
        particleTypeVector(particleInds) = 2;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(2,:),68,1);  %[kurze, lange]
        particle2List(sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 1 1], 68, 1);
    end
end

sumParticle = 0;    
%place horizontal particles
while sumParticle < numHorizontal
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(horizontalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = 1;
        particleTypeVector(particleInds) = 2;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(2,:),68,1);
        particle2List(numVertical+sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 1 1], 68, 1);
    end
end    
fprintf('%d particle 2 placed \n',numVertical+numHorizontal)
end
if i==3
%% place particle 3
fprintf('Start placing particle 3 (2x6): \n')
sumParticle = 0;
numVertical = floor(numParticle3/2);
% numVertical = 0;
numHorizontal = numParticle3-numVertical;
verticalInds = findParticleInds(6,2); % local index in stencil
horizontalInds = findParticleInds(2,6);

%place vertical particles
while sumParticle < numVertical
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(verticalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = -1;
        particleTypeVector(particleInds) = 3;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(3,:),12,1);  %[kurze, lange]
        particle3List(sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 -1 -1], 12, 1);
        
    end
end

sumParticle = 0;    
%place horizontal particles
while sumParticle < numHorizontal
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(horizontalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = -1;
        particleTypeVector(particleInds) = 3;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(3,:),12,1);
        particle3List(numVertical+sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([-1 -1 1 1], 12, 1);
    end
end    
fprintf('%d particle 3 placed \n',numVertical+numHorizontal)
end
if i == 4
%% place particle 4
fprintf('Start placing particle 4 (2x30): \n')
sumParticle = 0;
numVertical = floor(numParticle4/2);
numHorizontal = numParticle4-numVertical;
verticalInds = findParticleInds(30,2); % local index in stencil
horizontalInds = findParticleInds(2,30);

%place vertical particles
while sumParticle < numVertical
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(verticalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = -1;
        particleTypeVector(particleInds) = 4;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(4,:),60,1);  %[kurze, lange]
        particle4List(sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 -1 -1], 60, 1);
        
    end
end

sumParticle = 0;    
%place horizontal particles
while sumParticle < numHorizontal
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(horizontalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = -1;
        particleTypeVector(particleInds) = 4;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(4,:),60,1);
        particle4List(numVertical+sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([-1 -1 1 1], 60, 1);
    end
end    
fprintf('%d particle 4 placed \n',numVertical+numHorizontal)
end
if i == 5
%% place particle 5
fprintf('Start placing particle 5 (8x100): \n')
sumParticle = 0;
numVertical = floor(numParticle5/2);
numHorizontal = numParticle5-numVertical;
verticalInds = findParticleInds(100,8); % local index in stencil
horizontalInds = findParticleInds(8,100);

%place vertical particles
while sumParticle < numVertical
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(verticalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = -1;
        particleTypeVector(particleInds) = 5;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(5,:),800,1);  %[kurze, lange]
        particle5List(sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 -1 -1], 800, 1);
        
    end
end

sumParticle = 0;    
%place horizontal particles
while sumParticle < numHorizontal
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(horizontalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = -1;
        particleTypeVector(particleInds) = 5;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(5,:),800,1);
        particle5List(numVertical+sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([-1 -1 1 1], 800, 1);
    end
end    
fprintf('%d particle 5 placed \n',numVertical+numHorizontal)
end
if i == 6
%% place particle 6
fprintf('Start placing particle 6 (1x1): \n')
sumParticle = 0;
numVertical = floor(numParticle6/2);
numHorizontal = numParticle6-numVertical;
verticalInds = findParticleInds(1,1); % local index in stencil
horizontalInds = findParticleInds(1,1);

%place vertical particles
while sumParticle < numVertical
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(verticalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = 1;
        particleTypeVector(particleInds) = 6;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(2,:),1,1);  %[kurze, lange]
        particle6List(sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([1 1 -1 -1], 1, 1);
    end
end

sumParticle = 0;    
%place horizontal particles
while sumParticle < numHorizontal
    ind = randi(numSquares); % choose random position
    sten = stencil(g.NX,NZd,ind,stencilSize(i)); 
    particleInds = sten(horizontalInds); % get global inds
    if sum( bulkVector(particleInds) ~= 0 ) == 0
        sumParticle = sumParticle + 1;
        bulkVector(particleInds) = 1;
        bulkTypeVector(particleInds) = 1;
        particleTypeVector(particleInds) = 6;
        chargeInfo(particleInds,:) = repmat(chargeAtEdges(2,:),1,1);
        particle6List(numVertical+sumParticle,:) = particleInds;
        edgeChargeVector(g.CE0T(particleInds,:)) = repmat([-1 -1 1 1], 1, 1);
    end
end    
fprintf('%d particle 6 placed \n',numVertical+numHorizontal)
end

end

totalCharge = numParticle1*surfaceParticle1 + numParticle2*surfaceParticle2 - numParticle3*surfaceParticle3 - numParticle4*surfaceParticle4 - ...
    numParticle5*surfaceParticle5 - numParticle6*surfaceParticle6;
fprintf(' Total charge: %d \n',totalCharge)
end
