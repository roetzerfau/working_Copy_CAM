function [ globS ] = assembleGlobSsub( g , markE0Tint , hatSdiag , hatSoffdiag , eta , airVec , numAir )

K = numAir;
N = size(hatSdiag, 1);
globS = sparse(K*N, K*N);

for n = 1 : 4
    globS = globS + eta * ( kron( spdiags(markE0Tint(:,n), 0, K, K) , hatSdiag(:,:,n) ) ...
        - kron( g.markE0TE0Tsub{n}(airVec,airVec) , hatSoffdiag(:,:,n) ) );
end  %for

end  %function