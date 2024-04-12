 function edgeChargeVector = setEdgeChargeVector(g, bulkVector, edgeChargeVector)
    
     for i = 1 : g.numT
        if bulkVector(i) == 99
           possibleNeighbors = stencil( g.NX, g.NX, i, 1); 
           possibleNeighbors = possibleNeighbors(2:end); 
           for neigh = 1 : 4
              if ((bulkVector(possibleNeighbors(neigh)) == 0 ) ||(bulkVector(possibleNeighbors(neigh)) == 2 ))
                 if (neigh == 1)
                     edgeChargeVector(g.CE0T(i,1)) = 1;
                 elseif (neigh == 2)
                     edgeChargeVector(g.CE0T(i,4)) = 1;
                 elseif (neigh == 3)
                     edgeChargeVector(g.CE0T(i,3)) = 1;
                 elseif (neigh == 4)
                     edgeChargeVector(g.CE0T(i,2)) = 1;
                 end

              end
           end
        end
     end

 end