function [numEdgeTypesPOMparticles] =  countEdgeTypesPOMparticles(g, bulkVector, POMVector, POMParticleList, edgeChargeVector, reactiveSurfaceVector)

    numEdgeTypesPOMparticles = zeros(length(POMParticleList), 5);
    for i = 1 : length(POMParticleList)
        numPoreEdges = 0;
        numPOMedges = 0;
        numSolidStandardEdges = 0;
        numSolidReactiveEdges = 0;
        numSolidMemoryEdges = 0;
        for part = 1 : length(POMParticleList{i})
            particle = POMParticleList{i}(part);
            possibleNeighbors = stencil( g.NX, g.NX, particle, 1); 
            possibleNeighbors = possibleNeighbors(2:end); 
            
            for neigh = 1 : 4
                neighbor = possibleNeighbors(neigh);
                % only consider solid cells that have solid or POM
                % neighbor, not belonging to same POM particle
                if ( (bulkVector(neighbor) == 1) && ~ismember(neighbor, POMParticleList{i}) )
                    if (neigh == 1)
%                         particleEdgeInd = 1;
                        neighEdgeInd = 2;
                    elseif (neigh == 2) 
%                         particleEdgeInd = 4;
                        neighEdgeInd = 3;
                    elseif (neigh == 3)
%                         particleEdgeInd = 3;
                        neighEdgeInd = 4;
                    elseif (neigh == 4)
%                         particleEdgeInd = 2;
                        neighEdgeInd = 1;
                    end
                    
                    % POM - POM contact 
                    if (POMVector(neighbor) == 1)
                        numPOMedges = numPOMedges + 1;
                        
                    else
                        if ( edgeChargeVector(g.CE0T(neighbor,neighEdgeInd)) == 1)
                            numSolidMemoryEdges = numSolidMemoryEdges + 1;
                        elseif ( reactiveSurfaceVector(g.CE0T(neighbor,neighEdgeInd)) == 1)
                            numSolidReactiveEdges = numSolidReactiveEdges + 1;
                        else
                            numSolidStandardEdges = numSolidStandardEdges + 1;
                        end
                    end
                elseif ( bulkVector(neighbor) == 0 )
                    numPoreEdges = numPoreEdges + 1;
                end
            end
        end
        numEdgeTypesPOMparticles(i, :) = [numPoreEdges, numPOMedges, numSolidStandardEdges, ...
        numSolidReactiveEdges, numSolidMemoryEdges];
    end


end