function [numEdgeTypes] =  countEdgeTypes(g, bulkVector, POMVector, solidParticleList, edgeChargeVector, reactiveSurfaceVector, particleTypeVector)
    numSolidFluidEdges = 0;
    numSolidSolidStandardEdges = 0;
    numSolidSolidReactiveEdges = 0;
    numSolidPOMStandardEdges = 0;
    numSolidPOMReactiveEdges = 0;
    numSolidPOMMemoryEdges = 0;

    for i = 1 : length(solidParticleList)
        for part = 1 : length(solidParticleList{i})
            particle = solidParticleList{i}(part);
            possibleNeighbors = stencil( g.NX, g.NX, particle, 1); 
            possibleNeighbors = possibleNeighbors(2:end); 
            
            for neigh = 1 : 4
                neighbor = possibleNeighbors(neigh);
                % only consider solid cells that have solid or POM
                % neighbor, not belonging to same particle
                if ( (particleTypeVector(particle) ~= particleTypeVector(neighbor)) && (bulkVector(neighbor) == 1))
                    if (neigh == 1)
                        particleEdgeInd = 1;
                        neighEdgeInd = 2;
                    elseif (neigh == 2) 
                        particleEdgeInd = 4;
                        neighEdgeInd = 3;
                    elseif (neigh == 3)
                        particleEdgeInd = 3;
                        neighEdgeInd = 4;
                    elseif (neigh == 4)
                        particleEdgeInd = 2;
                        neighEdgeInd = 1;
                    end
                    
                    % solid - solid contact
                    if (POMVector(neighbor) == 0)
                        if ( (reactiveSurfaceVector(g.CE0T(particle,particleEdgeInd)) == 1 || ...
                              edgeChargeVector(g.CE0T(particle,particleEdgeInd)) == 1)  && ...
                             (reactiveSurfaceVector(g.CE0T(neighbor,neighEdgeInd)) == 1 || ...
                             edgeChargeVector(g.CE0T(neighbor,neighEdgeInd)) == 1))
                           numSolidSolidReactiveEdges = numSolidSolidReactiveEdges + 1; 
                        else
                            numSolidSolidStandardEdges = numSolidSolidStandardEdges + 1;
                        end
                    % solid - POM contact
                    else
                        if ( edgeChargeVector(g.CE0T(particle,particleEdgeInd)) == 1)
                            numSolidPOMMemoryEdges = numSolidPOMMemoryEdges + 1;
                        elseif( reactiveSurfaceVector(g.CE0T(particle,particleEdgeInd)) == 1)
                            numSolidPOMReactiveEdges = numSolidPOMReactiveEdges + 1;  
                        else
                            numSolidPOMStandardEdges = numSolidPOMStandardEdges + 1;
                        end
                    end
                elseif ( bulkVector(neighbor) == 0 )
                    numSolidFluidEdges = numSolidFluidEdges + 1;
                end
                
            end
        end
    end

% divide edges that were counted twice by two    
% numSolidSolidStandardEdges = numSolidSolidStandardEdges / 2;
% numSolidSolidReactiveEdges = numSolidSolidReactiveEdges / 2;
% numSolidSolidMemoryEdges = numSolidSolidMemoryEdges / 2;
numEdgeTypes = [numSolidFluidEdges, numSolidSolidStandardEdges, numSolidSolidReactiveEdges, ...
    numSolidPOMStandardEdges, numSolidPOMReactiveEdges, numSolidPOMMemoryEdges];

 
end
