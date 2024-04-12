function [ globKN ] = assembleGlobKNsub(  g , p , ord , markE0Tbdr , g_N , airVec , numAir )

global gPhi1D

K = numAir;
N = (p+1)^2;
[Q, W] = gaussQuadRule1D(ord);
Q2X = @(X,Y) g.deltaX * ones(K,1) * X + g.coordV0Tsub(airVec,1,1) * ones(size(X));
Q2Y = @(X,Y) g.BAsub(airVec) * X + g.DAsub(airVec) * Y + g.ACBDsub(airVec) * (X .* Y) + g.coordV0Tsub(airVec,1,2) * ones(size(X));

globKN = zeros(K, N);

for n = 1 : 4
    Length = g.lengthE0Tsub(:,n);
    Length = Length(airVec);
    switch n
        case 1, QX = Q;              QY = zeros(size(Q));
        case 2, QX = Q;              QY = ones(size(Q));
        case 3, QX = ones(size(Q));  QY = Q;
        case 4, QX = zeros(size(Q)); QY = Q;
    end  % siwtch
    gNAlgn = g_N( Q2X(QX.', QY.') , Q2Y(QX.', QY.') );
    Kkn = - markE0Tbdr(:,n) .* Length;
    for i = 1 : N
        globKN(:,i) = globKN(:,i) + Kkn .* ( gNAlgn * ( W .* gPhi1D(:,i,n).' ).' );
    end  % for i
end  % for n

globKN = reshape(globKN', K*N, 1);

end  % function