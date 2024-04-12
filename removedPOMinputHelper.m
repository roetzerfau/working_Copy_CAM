load('FinalConfig/config.1000.mat')
inputPOMtime = timeRemovedPOMparticles;
inputPOMparticles = removedPOMparticles;
inputPOMparticlesConc = removedPOMparticlesConc;
save('Input/inputPOM.mat','inputPOMtime','inputPOMparticles','inputPOMparticlesConc')
