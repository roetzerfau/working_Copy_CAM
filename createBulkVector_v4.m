function [bulkVector,bulkTypeVector,QuartzList, IlliteList, GoethitList,concAgent] = createBulkVector_v4( g , porosity , concAgent, IlliteCharge,QuartzCharge,GoethitCharge,NZd,QuartzRadius,GoethitLength, IlliteAmount, GoethitAmount,IlliteSize,agent)

rng('shuffle');
numSquares      = (g.NX)^2;
% borderInd       = [1 : g.NX ,(g.NX)^2-g.NX+1 : (g.NX)^2 ,  [1:g.NX-2]*g.NX+1 , [2:g.NX-1]*g.NX ];
% squarePool      = setdiff([1 : numSquares ] , borderInd , 'stable'); 
bulkVector      = zeros(numSquares, 1);
bulkTypeVector  = zeros(numSquares, 1);
numParticles    = floor(numSquares*(1-porosity)); 
% agent           = 1;


% numGoethit = floor(numParticles*0.05/GoethitLength);
numGoethit = floor(numParticles*GoethitAmount/GoethitLength);
numHorizontalGoethit = floor(numGoethit*0.5);
numVerticalGoethit = numGoethit-numHorizontalGoethit;

% numIllite  = floor(numParticles*0.20/6);
numIllite  = floor(numParticles*IlliteAmount/IlliteSize);
numHorizontalIllite = floor(numIllite*0.5);
numVerticalIllite = numIllite-numHorizontalIllite;

if QuartzRadius == 5
    load('circle_5.mat','circle');
elseif QuartzRadius == 10
    load('circle_10.mat','circle');
elseif QuartzRadius == 20
    load('circle_20.mat','circle');
elseif QuartzRadius == 50
    load('circle_50.mat','circle');
else
    assert('wrong Quartzradius')
end
numQuartz   = floor((numParticles - numGoethit*GoethitLength - numIllite*IlliteSize)/length(circle));

% numPosNeedles = floor(numParticles*numPos/3);
% numNegNeedles = floor(numParticles*numNeg/3);

% numVerticalPosNeedles   = floor(numPosNeedles*numVerticalPos);
% numVerticalNegNeedles   = floor(numNegNeedles*numVerticalNeg);
% numHorizontalPosNeedles = numPosNeedles - numVerticalPosNeedles; 
% numHorizontalNegNeedles = numNegNeedles - numVerticalNegNeedles;
% 
% needleList = zeros(numHorizontalNegNeedles+numHorizontalNegNeedles+numVerticalPosNeedles+numVerticalNegNeedles,3);
% horizontalNegNeedle = -1
% verticalNegNeedle   = -2
% horizontalPosNeedle = 1
% verticalPosNeedle   = 2
QuartzList = zeros(numQuartz,length(circle));
% IlliteList = zeros(numIllite,6);
IlliteList = zeros(numIllite,IlliteSize);
GoethitList = zeros(numGoethit,GoethitLength);

% needleCounter = 1;
% 0: 
QuartzRot  = zeros(numQuartz,1); 
IlliteRot  = zeros(numIllite,1);
GoethitRot = zeros(numGoethit,1);

%% place Quartz particles
sumQuartz = 0;
while sumQuartz<numQuartz
    ind = randi(numSquares);
    quartzStencil = stencil(g.NX,NZd,ind,g.NX);
%     quartzInds = setdiff(quartzStencil,quartzStencil([42 51 52 61]),'stable');
    quartzInds = quartzStencil(circle);

    if sum( bulkVector(quartzInds) ~= 0 ) == 0   
        bulkVector(quartzInds)  = 1;
        bulkTypeVector(quartzInds)   = QuartzCharge;

        sumQuartz = sumQuartz + 1;
        QuartzList(sumQuartz,:) = quartzInds;
        QuartzRot(sumQuartz)    = 1;

%         for bulk = quartzInds
%             concAgent(g.E0T(bulk,:)) = 0;
%         end
    end
end

fprintf('%d Quartz particles placed \n',sumQuartz)
%% place vertical Goethit particles
sumVerticalGoethit = 0;
while sumVerticalGoethit<numVerticalGoethit
    ind = randi(numSquares);
    goethitStencil = stencil(g.NX,NZd,ind,17);
    % 3  [1 2 5]
    % 5  [1 2 5 6 13]
    % 15 [1,2,5,6,13,14,25,26,41,42,61,62,85,86,113]
    if GoethitLength == 3
        goethitInds = [1 2 5];
    elseif GoethitLength == 5
        goethitInds = [1 2 5 6 13];
    elseif GoethitLength == 15
        goethitInds = [1,2,5,6,13,14,25,26,41,42,61,62,85,86,113];
    else 
        assert('wrong GoethitLength')
    end

    
    goethitInds    = goethitStencil(goethitInds); 
    if sum( bulkVector(goethitInds) ~= 0 ) == 0
        bulkVector(goethitInds)      = 1;
        bulkTypeVector(goethitInds)  = GoethitCharge;
        sumVerticalGoethit = sumVerticalGoethit + 1;
        GoethitList(sumVerticalGoethit,:) = goethitInds;
        GoethitRot(sumVerticalGoethit) = 1;

%         for bulk = goethitInds(1:end-2)
%             concAgent(g.E0T(bulk,:)) = [0 0 1 1];
%         end
%         concAgent(g.E0T(goethitInds(end-1),:)) = [1 0 1 1];
%         concAgent(g.E0T(goethitInds(end),:)) = [0 1 1 1];
    end
end
fprintf('%d vertical Goethit particles placed \n',sumVerticalGoethit)
%% place horizontal Goethit particles
sumHorizontalGoethit = 0;
while sumHorizontalGoethit<numHorizontalGoethit
    ind = randi(numSquares);
    goethitStencil = stencil(g.NX,NZd,ind,17);
    % 3  [1 3 4]
    % 5  [1 3 4 9 10]
    % 15 [1,3,4,9,10,19,20,33,34,51,52,73,74,99,100]

    if GoethitLength == 3
        goethitInds = [1 3 4];
    elseif GoethitLength == 5
        goethitInds = [1 3 4 9 10];
    elseif GoethitLength == 15
        goethitInds = [1,3,4,9,10,19,20,33,34,51,52,73,74,99,100];
    else 
        assert('wrong GoethitLength')
    end
    goethitInds    = goethitStencil(goethitInds); 
    if sum( bulkVector(goethitInds) ~= 0 ) == 0
        bulkVector(goethitInds)      = 1;
        bulkTypeVector(goethitInds)  = GoethitCharge;
        sumHorizontalGoethit = sumHorizontalGoethit + 1;
        GoethitList(sumVerticalGoethit+sumHorizontalGoethit,:) = goethitInds;
        GoethitRot(sumVerticalGoethit+sumHorizontalGoethit) = 2;

        
%         for bulk = goethitInds(1:end-2)
%             concAgent(g.E0T(bulk,:)) = [1 1 0 0];
%         end
%         concAgent(g.E0T(goethitInds(end-1),:)) = [1 1 0 1];
%         concAgent(g.E0T(goethitInds(end  ),:)) = [1 1 1 0];
    end
end
fprintf('%d horizontal Goethit particles placed \n',sumHorizontalGoethit)
%% place vertical Illite particles
sumVerticalIllite = 0;
while sumVerticalIllite<numVerticalIllite
    ind = randi(numSquares);
    illiteStencil = stencil(g.NX,NZd,ind,20);
    if IlliteSize == 6
        illiteInds    = illiteStencil([1 2 4 5 8 12]);
    elseif IlliteSize == 96
        illiteInds    = illiteStencil([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,53,55,56,57,58,59,60,61,63,64,65,66,67,68,69,77,79,80,81,82,83,84,85,89,90,91,92,93,105,107,108,109,110,111,112,119,120,121,137,139,140,141,142,153,173,175,176,213]);
    elseif IlliteSize == 24
        illiteInds    = illiteStencil([1,4,5,10,12,13,20,22,24,25,36,38,40,41,56,58,60,61,80,82,84,108,110,140]);
    end
    
    if sum( bulkVector(illiteInds) ~= 0 ) == 0
        bulkVector(illiteInds)      = 1;
        bulkTypeVector(illiteInds)  = IlliteCharge;
        sumVerticalIllite = sumVerticalIllite + 1;
        IlliteList(sumVerticalIllite,:) = illiteInds;

%         concAgent(g.E0T(illiteInds(1),:)) = [0 0 0 1];
%         concAgent(g.E0T(illiteInds(2),:)) = [1 0 0 1];
%         concAgent(g.E0T(illiteInds(3),:)) = [0 0 1 0];
%         concAgent(g.E0T(illiteInds(4),:)) = [0 1 0 1];
%         concAgent(g.E0T(illiteInds(5),:)) = [1 0 1 0];
%         concAgent(g.E0T(illiteInds(6),:)) = [0 1 1 0];
    end
end
fprintf('%d vertical Illite particles placed \n',sumVerticalIllite)
%% place horizontal Illite particles
sumHorizontalIllite = 0;
while sumHorizontalIllite<numHorizontalIllite
    ind = randi(numSquares);
    illiteStencil = stencil(g.NX,NZd,ind,20);
    if IlliteSize == 6
        illiteInds    = illiteStencil([1 3 4 5 11 12]); 
    elseif IlliteSize == 96
        illiteInds    = illiteStencil([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,67,68,69,70,71,72,73,75,76,77,78,79,80,81,82,93,94,95,96,97,101,103,104,105,106,107,108,123,124,125,133,135,136,137,138,157,169,171,172,209]);
    elseif IlliteSize == 24
        illiteInds    = illiteStencil([1,4,5,10,12,13,20,22,24,25,34,36,38,40,52,54,56,58,76,78,80,104,106,136]);
    end

    if sum( bulkVector(illiteInds) ~= 0 ) == 0
        bulkVector(illiteInds)      = 1;
        bulkTypeVector(illiteInds)  = IlliteCharge;
        sumHorizontalIllite = sumHorizontalIllite + 1;
        IlliteList(sumVerticalIllite+sumHorizontalIllite,:) = illiteInds;

%         concAgent(g.E0T(illiteInds(1),:)) = [1 0 0 0];
%         concAgent(g.E0T(illiteInds(2),:)) = [1 0 0 1];
%         concAgent(g.E0T(illiteInds(3),:)) = [1 0 1 0];
%         concAgent(g.E0T(illiteInds(4),:)) = [0 1 0 0];
%         concAgent(g.E0T(illiteInds(5),:)) = [0 1 0 1];
%         concAgent(g.E0T(illiteInds(6),:)) = [0 1 1 0];
    end
end
fprintf('%d horizontal Illite particles placed \n',sumHorizontalIllite)
end