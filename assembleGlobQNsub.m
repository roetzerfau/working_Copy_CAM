function [ globQN ] = assembleGlobQNsub( g , markE0Tbdr , hatQdiag , airVector , numAir )

K = numAir;
N = size(hatQdiag, 1);
globQN = cell(2, 1);
globQN{1} = sparse(K*N, K*N);
globQN{2} = sparse(K*N, K*N);

for n = 1 : 4
    NuLength1 = g.NuLengthE0Tsub(:,n,1);
    NuLength1 = NuLength1(airVector);
    NuLength2 = g.NuLengthE0Tsub(:,n,2);
    NuLength2 = NuLength2(airVector);
    globQN{1} = globQN{1} + kron( spdiags(markE0Tbdr(:,n) .* NuLength1, ...
        0, K, K) , hatQdiag(:,:,n) );
    globQN{2} = globQN{2} + kron( spdiags(markE0Tbdr(:,n) .* NuLength2, ...
        0, K, K) , hatQdiag(:,:,n) );
end  % for

end  % function