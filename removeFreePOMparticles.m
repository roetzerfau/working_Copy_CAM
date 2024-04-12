 function [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, POMParticleList,...
     removedPOMparticles, removedPOMparticlesConc, timeRemovedPOMparticles, totalPOMoutputConc] = ...
     removeFreePOMparticles(g, bulkVector, bulkTypeVector, edgeChargeVector, reactiveSurfaceVector,...
     POMVector, POMconcVector, POMageVector, POMParticleList, removePOMthreshold, freePOMparticles,...
     removedPOMparticles, removedPOMparticlesConc, timeRemovedPOMparticles, totalPOMoutputConc, k)
    
    if(~isempty(freePOMparticles))
        indRemovePOMparticles = freePOMparticles(find(freePOMparticles(:,2)>removePOMthreshold),1);

        for i = 1 : length(indRemovePOMparticles)
            removedPOMparticles = [removedPOMparticles, POMParticleList{indRemovePOMparticles(i)}];
            timeRemovedPOMparticles = [timeRemovedPOMparticles k];
            removedPOMparticlesConc = [removedPOMparticlesConc, POMconcVector(POMParticleList{indRemovePOMparticles(i)})];
            POMVector(POMParticleList{indRemovePOMparticles(i)}) = 0;
            bulkVector(POMParticleList{indRemovePOMparticles(i)}) = 0;
            bulkTypeVector(POMParticleList{indRemovePOMparticles(i)}) = 0;
            totalPOMoutputConc = totalPOMoutputConc + sum(POMconcVector(POMParticleList{indRemovePOMparticles(i)}));
            POMconcVector(POMParticleList{indRemovePOMparticles(i)}) = 0;
            POMageVector(POMParticleList{indRemovePOMparticles(i)}) = 0;
        end

        for i = length(indRemovePOMparticles): -1 : 1
            POMParticleList(indRemovePOMparticles(i)) = [];
        end
    end

 end