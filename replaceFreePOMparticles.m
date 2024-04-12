 function [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList] = ...
     replaceFreePOMparticles(g, bulkVector, bulkTypeVector, edgeChargeVector, ...
     reactiveSurfaceVector, POMVector, POMconcVector, POMageVector, POMParticleList)
    
    POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);
    
    indFreePOMparticles = [];
    
    for i = 1 : length(POMParticleList)
        if sum(edgeChargeVector(POMsolidEdgeList{i})) + sum(reactiveSurfaceVector(POMsolidEdgeList{i})) == 0
           indFreePOMparticles = [indFreePOMparticles i];
        end
    end
 
    indFreePOMparticles
    
    movingPOMparticles = cell(1, length(indFreePOMparticles));
    
    for i = 1 : length(indFreePOMparticles)
        movingPOMparticles{i} = POMParticleList{indFreePOMparticles(i)};
        POMVector(POMParticleList{indFreePOMparticles(i)}) = 0;
        bulkVector(POMParticleList{indFreePOMparticles(i)}) = 0;
        bulkTypeVector(POMParticleList{indFreePOMparticles(i)}) = 0;
        POMconcVector(POMParticleList{indFreePOMparticles(i)}) = 0;
        POMageVector(POMParticleList{indFreePOMparticles(i)}) = 0;
    end
    
    for i = length(indFreePOMparticles): -1 : 1
        POMParticleList(indFreePOMparticles(i)) = [];

    end
    
    for i = 1 : length(indFreePOMparticles)
         [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList] = ...
            placePOMparticleRandomly(g, bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList, ...
            {movingPOMparticles{i}}, length(movingPOMparticles{i}));
    end

 end