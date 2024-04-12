 function [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, totalPOMinputConc] = ...
     placePOMparticleRandomly(g, bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, ...
     randomPOMparticles, randomPOMparticlesSizes, totalPOMinputConc)
    
    distFromSolidThreshold = 0;
    startConcPOM = 1;
 
    n = g.NX;
 
    geo = reshape(bulkVector, [n n]);
    geo = repmat(geo, 3);
    
    [distMatrix,~] = bwdist(geo,'cityblock');
    geo = geo(n + 1 : 2*n, n + 1 : 2*n);
    distMatrix = distMatrix(n + 1 : 2*n, n + 1 : 2*n);
    
    geo = reshape(geo, [n*n 1]);
    distMatrix = reshape(distMatrix, [n*n 1]);
    
    
    numPOMparticles = length(randomPOMparticles);
    sizePOMparticles = zeros(numPOMparticles,1);
    for i = 1 : numPOMparticles
        sizePOMparticles(i) = length(randomPOMparticles{i});
    end
    % if particles of certain size should be placed for initial state
%     cand = find(sizePOMparticles == 10);
%     candInd = randi(length(cand));
%     indPOMparticle = cand(candInd);
    indPOMparticle = randi(numPOMparticles);
    
    stencilPOMparticle = stencil(n, n, randomPOMparticles{indPOMparticle}(1), max(2*(randomPOMparticlesSizes(indPOMparticle)-1),1));
    indHelperPOMparticle = ismember(stencilPOMparticle, randomPOMparticles{indPOMparticle});
    
    newPositionFound = 0;
    while ~newPositionFound
        % if particles of certain size should be placed for initial state
%         potentialAims = find(distMatrix > max(randomPOMparticlesSizes(indPOMparticle) - 4,1));
%         helper = randi(length(potentialAims));
%         aimInd = potentialAims(helper);
%         stencilNewPosition = stencil(n, n, aimInd, max(2*(randomPOMparticlesSizes(indPOMparticle)-1),1));
        stencilNewPosition = stencil(n, n, randi(n*n), max(2*(randomPOMparticlesSizes(indPOMparticle)-1),1));
        globalIndNewPOMparticle = stencilNewPosition(indHelperPOMparticle);
        if min(min(distMatrix(globalIndNewPOMparticle))) > distFromSolidThreshold
            newPositionFound = 1;
        end
    end
    bulkVector(globalIndNewPOMparticle) = 1;
    POMVector(globalIndNewPOMparticle) = 1;
    bulkTypeVector(globalIndNewPOMparticle) = -1;
    POMconcVector(globalIndNewPOMparticle) = 1;
    POMageVector(globalIndNewPOMparticle) = 1;
    POMParticleList{length(POMParticleList) + 1} = globalIndNewPOMparticle;
    totalPOMinputConc = totalPOMinputConc + length(globalIndNewPOMparticle);
 end
