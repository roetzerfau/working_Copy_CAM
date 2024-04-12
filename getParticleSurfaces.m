function [particleSurface] = getParticleSurfaces(g, particleList, particleTypeVector)

   particleSurface = zeros(length(particleList),1);
   
   for i = 1 : length(particleList)
       for j = 1 : length(particleList{i})
           solidInd = particleList{i}(j);
           possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
           possibleNeighbors = possibleNeighbors(2:end);
           for neigh = 1 : 4
                if (neigh == 1)
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurface(i) = particleSurface(i) + 1;
                    end
                elseif (neigh == 2)
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurface(i) = particleSurface(i) + 1;
                    end
                elseif (neigh == 3)
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurface(i) = particleSurface(i) + 1;
                    end
                elseif (neigh == 4)
                    if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                        particleSurface(i) = particleSurface(i) + 1;
                    end
                end
            end
       end
   end
end

    