function [ globR ] = assembleGlobRsub( g , markE0Tint , hatRdiag , hatRoffdiag , ...
    K11DG , K12DG , K21DG , K22DG , airVec , numAir )

[~, N] = size(K11DG);
K = numAir;
globR = cell(2, 1);
globR{1} = sparse(K*N, K*N);
globR{2} = sparse(K*N, K*N);

for n = 1 : 4
    NuLength1 = g.NuLengthE0Tsub(:,n,1);
    NuLength1 = NuLength1(airVec);
    NuLength2 = g.NuLengthE0Tsub(:,n,2);
    NuLength2 = NuLength2(airVec);
    for i = 1 : N
        globR{1} = globR{1} + kron( bsxfun(@times, g.markE0TE0Tsub{n}(airVec,airVec), (0.5 * ...
            NuLength1 .* K11DG(:,i) + 0.5 * NuLength2 .* ...
            K21DG(:,i)).' ) , hatRoffdiag(:,:,i,n) ) ...
            + kron( spdiags(0.5 * markE0Tint(:,n) .* NuLength1 .* K11DG(:,i) ...
            + 0.5 * markE0Tint(:,n) .* NuLength2 .* K21DG(:,i), 0, K, K) , ...
            hatRdiag(:,:,i,n) );
        globR{2} = globR{2} + kron( bsxfun(@times, g.markE0TE0Tsub{n}(airVec,airVec), (0.5 * ...
            NuLength1 .* K12DG(:,i) + 0.5 * NuLength2 .* ...
            K22DG(:,i)).' ) , hatRoffdiag(:,:,i,n) ) ...
            + kron( spdiags(0.5 * markE0Tint(:,n) .* NuLength1 .* K12DG(:,i) ...
            + 0.5 * markE0Tint(:,n) .* NuLength2 .* K22DG(:,i), 0, K, K) , ...
            hatRdiag(:,:,i,n) );
    end  % for i
end  % for n

end  % function