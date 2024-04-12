function [ conc ] = computeConcAir( g , uDG , ord )

global gPhi2D

[~, ~, W] = gaussQuadRule2D(ord);

conc = uDG * gPhi2D' * W' .* g.areaT;

end  % function