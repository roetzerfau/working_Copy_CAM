 function  POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList)
    
    POMsolidEdgeList = cell(1,length(POMParticleList));
    for i = 1 : length(POMParticleList)
        for part = 1 : length(POMParticleList{i})
            possibleNeighbors = stencil( g.NX, g.NX, POMParticleList{i}(part), 1); 
            possibleNeighbors = possibleNeighbors(2:end); 
            
            for neigh = 1 : 4
                if ((neigh == 1) && (bulkVector(possibleNeighbors(neigh)) == 1) &&  (POMVector(possibleNeighbors(neigh)) == 0))
                    POMsolidEdgeList{i} = [POMsolidEdgeList{i}; g.CE0T(possibleNeighbors(neigh),2)];
                elseif ((neigh == 2) && (bulkVector(possibleNeighbors(neigh)) == 1) &&  (POMVector(possibleNeighbors(neigh)) == 0))
                    POMsolidEdgeList{i} = [POMsolidEdgeList{i}; g.CE0T(possibleNeighbors(neigh),3)];
                elseif ((neigh == 3) && (bulkVector(possibleNeighbors(neigh)) == 1) &&  (POMVector(possibleNeighbors(neigh)) == 0))
                    POMsolidEdgeList{i} = [POMsolidEdgeList{i}; g.CE0T(possibleNeighbors(neigh),4)];
                elseif ((neigh == 4) && (bulkVector(possibleNeighbors(neigh)) == 1) &&  (POMVector(possibleNeighbors(neigh)) == 0))
                    POMsolidEdgeList{i} = [POMsolidEdgeList{i}; g.CE0T(possibleNeighbors(neigh),1)];
                end
            end
        end
    end

 end