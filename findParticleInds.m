function particleInds = findParticleInds(height,width)
% finds the inds of a stencil corresponding to the particle

% [ stencil ] = stencil( NX , NZd , candidates , layers )
% domain is 250x250 (can be too small, if particle is large enough)

if(mod(height,2)==1) % height is odd
    top = 125+(height-1)/2;
    bottom = 125-(height-1)/2;    
else
    top = 125+height/2;
    bottom = 125-(height/2-1);
end

if(mod(width,2)==1)
    left = 125-(width-1)/2;
    right = 125+(width-1)/2;
else
    left = 125-(width/2);
    right = 125+(width/2-1);
end

globalInds = [];
for row = bottom:top
for col = left:right
     globalInds = [globalInds, convertIndex([row,col],250)];
end
end
% sten contains stencil of size 100 for 250x250 geo around the center
sten = stencil(250,250,convertIndex([125,125],250),100);
% find local inds
particleInds = [];
for globalInd = globalInds
    particleInds = [particleInds,find(sten==globalInd)];
end

% particleInds contains indices of particle in stencil
particleInds = sort(particleInds);

end

function ind = convertIndex(coordinates,NX)
x = coordinates(1);
y = coordinates(2);
ind = (x-1)*NX+y;
end