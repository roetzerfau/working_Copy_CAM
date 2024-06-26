function [ globH ] = assembleGlobHsub( g , hatHc , hatHx , hatHy , airVec , numAir )

K = numAir; N = size(hatHc, 1);
globH = cell(2, 1);
globH{1} = sparse(K*N, K*N);
globH{2} = sparse(K*N, K*N);

globH{1} =   kron( spdiags(g.BAsub(airVec)  , 0, K, K) , hatHc(:,:,2) ) ...
           + kron( spdiags(g.ACBDsub(airVec), 0, K, K) , hatHy(:,:,2) ) ...
           - kron( spdiags(g.DAsub(airVec)  , 0, K, K) , hatHc(:,:,1) ) ...
           - kron( spdiags(g.ACBDsub(airVec), 0, K, K) , hatHx(:,:,1) );
globH{2} = - kron( spdiags(g.deltaX*ones(K,1), 0, K, K) , hatHc(:,:,2) );

end  % function