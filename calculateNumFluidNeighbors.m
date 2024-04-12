 function numFluidNeighVector = calculateNumFluidNeighbors(g, bulkVector, edgeChargeVector, POMVector, flag)
 % flag == 1: set numFluidNeigh to zero if POM particle has no attractive edge neighbors   
 numFluidNeighVector = zeros(length(bulkVector), 1);
 
     for i = 1 : g.numT
        numFluidNeigh = 0;
        if POMVector(i) == 1
           possibleNeighbors = stencil( g.NX, g.NX, i, 1); 
           possibleNeighbors = possibleNeighbors(2:end); 
           for neigh = 1 : 4
              if (bulkVector(possibleNeighbors(neigh)) == 0 )
                  numFluidNeigh = numFluidNeigh + 1;  
              end
           end
           numFluidNeighVector(i) = numFluidNeigh;
        end
     end
 
 end