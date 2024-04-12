for rep = 1 : 50
    load(strcat('FinalConfig/config.', int2str(rep), '.mat'))
    indFreePOMparticles = [];

    POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);
    
    for i = 1 : length(POMParticleList)
        if sum(edgeChargeVector(POMsolidEdgeList{i})) + sum(reactiveSurfaceVector(POMsolidEdgeList{i})) == 0
           indFreePOMparticles = [indFreePOMparticles i];
        end
    end

    indFreePOMparticles
    numFreePOMparticles = length(indFreePOMparticles);

end