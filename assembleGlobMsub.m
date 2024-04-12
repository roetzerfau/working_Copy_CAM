function [ globM ] = assembleGlobMsub( g , hatMc , hatMx , airVec , numAir )

globM = kron( spdiags(g.deltaX * g.DAsub(airVec), 0, numAir, numAir) , hatMc ) ...
    + kron( spdiags(g.deltaX * g.ACBDsub(airVec), 0, numAir, numAir) , hatMx );

end  % function