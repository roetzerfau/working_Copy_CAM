function MainSlim

%% Compiling the C++ Components, necessary to determine the particle size distribution
% mex particleSizeDistribution.cpp;

%% Definition of parameters                %% Can be changed
bilder = 1;                                % 0: nur Bild von Anfangs- und Endzustand, 1: Bild zu jedem Schritt, 2: Bild zu ersten 10 Zuständen und Endzustand, 3: Bild nach allen 100 Schritten

attraction_type = 1; % 1: uniform, 2: charges

% Number of Time Steps
numOuterIt  = 100;     

% Parameters
N = 100;
porosity = 0.5;

%% Creating domain for Simulation
[g, bulkVector, bulkTypeVector, reactiveSurfaceVector, particleTypeVector,...
    solidParticleList] = initializeDomainRandom(N, porosity);
  
NZd = g.NX; 

uDG(find(bulkVector == 0),1) = 0;
uDG(find(bulkVector == 1),1) = 1;
uDG(find((bulkVector == 1) & (particleTypeVector == 0)),1) = 2;
    
uLagr       = projectDG2LagrangeSub( uDG );
visualizeDataSub(g, uLagr, 'u', 'solu', 0);

fileID = fopen( 'Move_bulk_log_file' , 'w' );
% fileID_1 = fopen( 'Verteilungen', 'w' );
%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed
% printInfoSlim(0,bulkVector);

for k = 1 : numOuterIt
%% Doing the movement of Bulk
T_start = tic;
% move solid particles

for solidParticle = 1 : length( solidParticleList )
    particleSize = length( solidParticleList{ solidParticle } ); 

    bigParticleStencilLayers_individual = 5;
%     bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
%     solidParticle
    [bulkVector,bulkTypeVector, particleTypeVector,reactiveSurfaceVector,...
        solidParticleList{ solidParticle },~] = moveParticlesSlim( particleSize, bigParticleStencilLayers_individual,...
        g, bulkVector, bulkTypeVector, particleTypeVector,reactiveSurfaceVector, ...
        NZd , fileID ,solidParticleList{ solidParticle },4,0, attraction_type);  

end

% tic

% %% move big particles, typeflag = 4
% % T_start = tic;
% bigJumping = 1;
% if bigJumping == 1
% [particleList, particleContent] = particleInfoSlim(bulkVector, solidParticleList);
% for particle = 1 : length( particleList )
%     particleSize = length( particleList{ particle } ); 
%     if(size(particleContent{particle},1)<2)% kein Verbund
%         continue
%     end
%     test_ind = particleList{particle}(1);
%     
% %     bigParticleStencilLayers_individual = round(bigParticleStencilLayersConstant*1/(particleSize)^0.5);
%     bigParticleStencilLayers_individual = min(5, ceil(20/(particleSize)^0.5));
% % bigParticleStencilLayers_individual = 1;
%     if particleSize > 20000
%         bigParticleStencilLayers_individual = 0;
%     end
%     if bigParticleStencilLayers_individual == 0
%         continue
%     end
%     
%     [bulkVector,bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
%         concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, changedList,~] = ...
%         moveParticles( particleSize, bigParticleStencilLayers_individual, g, bulkVector, bulkTypeVector, ... 
%         particleTypeVector, POMVector, POMconcVector, POMageVector, concAgent, ...
%         concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, NZd , ...
%         fileID ,particleList{ particle },sumAgent,4,0, attraction_type);  
%   
% end
% end
% 
% fprintf('Time for jumping of aggregates: %d ', toc(T_start))
% % size(find(particleTypeVector))
% % size(find(bulkVector))
%%
T2 = tic;

    uDG(find(bulkVector == 0),1) = 0;
    uDG(find(bulkVector == 1),1) = 1;
    uDG(find((bulkVector == 1) & (particleTypeVector == 0)),1) = 2;
    
    if bilder == 0 && k == numOuterIt
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k);
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k);
    elseif bilder == 1 && (k <= 5 || mod(k,100) == 0 || k == numOuterIt)
% elseif bilder == 1 
    uLagr       = projectDG2LagrangeSub( uDG );
    visualizeDataSub(g, uLagr, 'u', 'solu', k);
    visualizeDataEdges(g, reactiveSurfaceVector, 'reactiveEdges', 'reactiveSurfaceVector', k, 2);
    elseif bilder == 2 && (k <= 10 || k == numOuterIt)	
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k); 
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k);  
    elseif bilder == 3 && (mod(k,100) == 0 || k == numOuterIt)	
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k); 
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k); 
    end

    %bacLagr     = projectDG2LagrangeSub( uDGbac );
    %visualizeDataSub(g, bacLagr, 'bac', 'solBac', k);
    %biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
    %visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', k);
    
%     particleDistribution = particleSizeDistribution(bulkVector.');
    


% printInfo(k,bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List, particle6List, particle7List, particle8List, ...
%     particle9List, particle10List, particle11List, particle12List, particle13List, particle14List, particle15List, particle16List, particle17List, ...
%     particle18List, particle19List, particle20List)

%     numberOfSolidParticles = sum(particleDistribution(:,1) .* particleDistribution(:,2));

%     if sum(bioMvector .* bulkVector ~= 0)
%         error('NEEEEEEEEEEEIN!')
%     end
%     particleDistribution

% printInfoTUM(k,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
%     numFreePOMparticles, numEdgeTypes, totalPOMinputConc, sumExcessPOM, POMagentInput);
% %     if(k < 25 || mod(k,25) == 0 || k == numOuterIt)
%     if(k < 5 || mod(k,50) == 0 || k == numOuterIt)
%         particleListHelper = particleList;
%         particleList = solidParticleList;
%         fileName    = ['FinalConfig/config','.', num2str(k),'.mat']; 
%         save(fileName,'g','bulkVector','bulkTypeVector','POMconcVector', 'concPOMAgent','edgeChargeVector','POMagentAge',...
%             'POMVector', 'POMageVector', 'POMParticleList', 'particleList', 'reactiveSurfaceVector', 'particleTypeVector')        
%         particleList = particleListHelper;
%     end
    
end  % for k
fclose( fileID );
% fclose( fileID_1);
end
