function particleList = createSolidParticleList (particleTypeVector) 

    numSolidParticles = max(max(particleTypeVector));
    particleList = cell(1,numSolidParticles);
    
    for i = 1 : numSolidParticles
        indSolidParticle = find(particleTypeVector == i);
        particleList{i} = indSolidParticle';
    end
end
