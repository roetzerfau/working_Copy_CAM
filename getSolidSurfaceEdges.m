function [particleSurfaceEdgeList] = getSolidSurfaceEdges(g, solidParticleList, particleTypeVector)
   particleSurfaceEdgeList = cell(length(solidParticleList),1);
   
   for i = 1 : length(solidParticleList)
       for j = 1 : length(solidParticleList{i})
           solidInd = solidParticleList{i}(j);
           possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
           possibleNeighbors = possibleNeighbors(2:end);
           for neigh = 1 : 4
                if (neigh == 1)
                    edgeDirection = 1;
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurfaceEdgeList{i} = [particleSurfaceEdgeList{i} g.CE0T(solidInd, edgeDirection)]; 
                    end
                elseif (neigh == 2)
                    edgeDirection = 4;
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurfaceEdgeList{i} = [particleSurfaceEdgeList{i} g.CE0T(solidInd, edgeDirection)];
                    end
                elseif (neigh == 3)
                    edgeDirection = 3;
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurfaceEdgeList{i} = [particleSurfaceEdgeList{i} g.CE0T(solidInd, edgeDirection)];
                    end
                elseif (neigh == 4)
                    edgeDirection = 2;
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurfaceEdgeList{i} = [particleSurfaceEdgeList{i} g.CE0T(solidInd, edgeDirection)];
                    end
                end
            end
       end
   end
end

    