function area = contact_area(bulkVector, particle1List, particle2List, particle3List, particle4List, particle5List,particle6List, particleList)
    surface = particle_surface(bulkVector,particleList);
    area = (size(particle1List,1)*36 + size(particle2List,1)*72+ size(particle3List,1)*16 + size(particle4List,1)*64 + ...
        size(particle5List,1)*216 + size(particle6List,1)*4 - surface) / 2;
end