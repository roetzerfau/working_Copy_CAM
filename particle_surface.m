function surf = particle_surface(bulkVector, particleList)
surf = 0;
NX = sqrt(length(bulkVector));
for particle = 1 : length( particleList )
    for bulk = particleList{ particle }
        sten = stencil( NX , NX , bulk , 1 );
        surf = surf + 4 - sum(bulkVector( sten(2:end) ) );
    end
end
end