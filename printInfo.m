function printInfo(k,bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List,particle6List)

particleListHelper = {particle1List,particle2List,particle3List,particle4List,particle5List,particle6List};
numParticles = zeros(1,6);
for i = 1 : 6
    numParticles(i) = size(particleListHelper{i},1); 
end
indParticles = find(numParticles > 0);
if size(indParticles,2) == 1
    indParticles(2) = indParticles(1);
end

particleDistribution = particleSizeDistribution(bulkVector.');
writematrix(particleDistribution,['Distribution/particleDistribution','.',num2str(k),'.txt'],'Delimiter','tab')
[particleList, particleContent] = particleInfoNew(bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List,particle6List);
fileName    = ['Distribution/Verteilung','.', num2str(k)];
fileID_k = fopen(fileName,'w');
fprintf( fileID_k , 'Kenngroessen für alle Partikel:\n');
fprintf( fileID_k , 'Anzahl der Teilchen : %d\n',length(particleList));
fprintf( fileID_k , 'mittlere Partikelgr.: %.3f\n', (particleDistribution(:,1)'*particleDistribution(:,2))/sum(particleDistribution(:,2))  );
fprintf( fileID_k , 'mean hydrodynamic diameter: %.3f\n', (2*sqrt(particleDistribution(:,1)/pi)'*particleDistribution(:,2))/sum(particleDistribution(:,2))  );
fprintf( fileID_k , 'Anzahl der Teilchen== Teilchen 1: %d\n',sum(cellfun('length',particleList)== size(particleListHelper{indParticles(1)},2)));
fprintf( fileID_k , 'Anzahl der Teilchen== Teilchen 2: %d\n',sum(cellfun('length',particleList)== size(particleListHelper{indParticles(2)},2)));
fprintf( fileID_k , 'Volumen: %d\n',sum(bulkVector));
fprintf( fileID_k , 'Oberfläche Gesamt: %d\n',particle_surface(bulkVector,particleList));
fprintf( fileID_k , 'Oberfläche^2/Volumen: %.3f\n',particle_surface(bulkVector,particleList)^2/sum(bulkVector));
kontaktflaeche = contact_area(bulkVector, particle1List, particle2List, particle3List, particle4List, particle5List, particle6List,particleList);
fprintf( fileID_k , 'Kontaktfläche Gesamt: %d\n',kontaktflaeche);
fprintf( fileID_k , 'Kontaktfläche durch Volumen: %d\n',kontaktflaeche/sum(bulkVector));
[particleList_loecher, ~] = particleInfoNew(~bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List,particle6List,1);
fluid_zshkmp = length(particleList_loecher);
[particleList_loecher, ~] = particleInfoNew(bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List,particle6List);
solid_zshkmp = length( particleList_loecher );
loecher = solid_zshkmp - fluid_zshkmp + 1;
fprintf( fileID_k , 'Loecher: %d\n', loecher);
fprintf( fileID_k , 'Volumengemittelte mittlere Partikelgr.: %.3f\n\n',1/sum(bulkVector)*((particleDistribution(:,1)').^2*particleDistribution(:,2)));

fprintf( fileID_k , 'Kenngroessen für Partikelverbuende:\n');
tmp = 0;
bulkVector_zsh = bulkVector;
particle_del_list = [];
for i = 1:length(particleList)
    if size(particleContent{i},1)>1
        tmp = tmp + 1;        
    else
        particle_del_list = [particle_del_list, i];
    end
end
for i = 1:length(particle_del_list)
   bulkVector_zsh(particleList{particle_del_list(i)}) = 0;
end
zsh_particle_ind = setdiff([1:length(particleList)],particle_del_list);
particleList_zsh = cell(1,length(zsh_particle_ind));

for i = 1:length(zsh_particle_ind)
    particleList_zsh{i}=particleList{zsh_particle_ind(i)};
end
   
   
fprintf( fileID_k , 'Anzahl der Teilchen : %d\n',tmp);
particleDistribution = particleSizeDistribution(bulkVector_zsh.');
fprintf( fileID_k , 'mittlere Partikelgr.: %.3f\n', (particleDistribution(:,1)'*particleDistribution(:,2))/sum(particleDistribution(:,2))  );
fprintf( fileID_k , 'mean hydrodynamic diameter: %.3f\n', (2*sqrt(particleDistribution(:,1)/pi)'*particleDistribution(:,2))/sum(particleDistribution(:,2))  );
fprintf( fileID_k , 'Anzahl der Teilchen>= Quartz: %d\n',sum(cellfun('length',particleList_zsh)>=0));
fprintf( fileID_k , 'Volumen: %d\n',sum(bulkVector_zsh));
fprintf( fileID_k , 'Oberfläche Gesamt: %d\n',particle_surface(bulkVector_zsh,particleList_zsh));
fprintf( fileID_k , 'Oberfläche^2/Volumen: %.3f\n',particle_surface(bulkVector_zsh,particleList_zsh)^2/sum(bulkVector_zsh)); 
[particleList_loecher, ~] = particleInfoNew(~bulkVector_zsh,particle1List,particle2List,particle3List,particle4List,particle5List,particle6List,1);
fluid_zshkmp = length(particleList_loecher);
[particleList_loecher, ~] = particleInfoNew(bulkVector_zsh,particle1List,particle2List,particle3List,particle4List,particle5List,particle6List);

solid_zshkmp = length( particleList_loecher );
loecher = solid_zshkmp - fluid_zshkmp + 1;
fprintf( fileID_k , 'Loecher: %d\n', loecher);
fprintf( fileID_k , 'Volumengemittelte mittlere Partikelgr.: %.3f\n',1/sum(bulkVector_zsh)*((particleDistribution(:,1)').^2*particleDistribution(:,2)));

fclose(fileID_k);

end


