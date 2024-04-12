function [ hatInt ] = computeHatInt( p , ord )

global gPhi2D

numBases = (p+1)^2;
[~, ~, W] = gaussQuadRule2D(ord);
hatInt = zeros(numBases,numBases,numBases);

for i = 1 : numBases
    for j = 1 : i
        for k = 1 : j
            hatInt(i,j,k) = W * ( gPhi2D(:,i) .* gPhi2D(:,j) .* gPhi2D(:,k) );
            hatInt(i,k,j) = hatInt(i,j,k);
            hatInt(j,i,k) = hatInt(i,j,k);
            hatInt(j,k,i) = hatInt(i,j,k);
            hatInt(k,i,j) = hatInt(i,j,k);
            hatInt(k,j,i) = hatInt(i,j,k);
        end  % for k
    end  % for j
end  % for i

end  % function