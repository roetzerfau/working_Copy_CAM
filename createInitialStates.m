numSimulations = 10;

%% Definition of parameters                %% Can be changed
NX          = 500;                          % Number of horizontal Elements
NZd         = NX;                          % Number of vertical Elements
porosity    = 0.9;                         % Porosity

% nicht relevant
chargeAtEdges = [-1,1;... % particle1 1x17  (Goethit)
                 -1,1;... % particle2 2x34  (Goethit)
                 -1,1;... % particle3 2x6   (Illite)
                 -1,1;... % particle4 2x30  (Illite)
                 -1,1;... % particle5 8x100 (Illite)
                 -1,1;... % particle6 1x12 (Illite)
                 -1,1;... % particle7 3x4 (Illite)
                 -1,1;... % particle8 1x36 (Illite)
                 -1,1;... % particle9 1x30 (Illite)
                 -1,1;... % particle10 2x15 (Illite)
                 -1,1;... % particle11 3x10 (Illite)
                 -1,1;... % particle12 5x6 (Illite)
                 -1,1;... % particle13 1x60 (Illite)
                 -1,1;... % particle14 3x20 (Illite)
                 -1,1;... % particle15 4x15 (Illite)
                 -1,1;... % particle16 5x12 (Illite)
                 -1,1;... % particle17 1x48 (Illite)
                 -1,1;... % particle18 2x24 (Illite)
                 -1,1;... % particle19 3x16 (Illite)
                 -1,1];  % particle20 6x8 (Illite)     
                 
amountOfParticles = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]; %Amount of particle 1, particle 2, ... (sum should be 1)

%% Fixed Parameters (of general framework)
NZu         = 0;
width       = NX;
intBound    = NZd * ones(NX+1,1);
upBound     = NZd * ones(NX+1,1);

%% Creating domain for Simulation
g = createDomainFolded(width, NX, NZd, NZu, intBound, upBound);

NX = g.NX;

for i=1:numSimulations
    %% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed
    concAgent       = 0*ones( g.numE , 1 );    % Initial Concentration of Agent
    edgeChargeVector = 0*ones( g.numCE , 1);    % Initial charge Vector
    
    [bulkVector,bulkTypeVector,chargeInfo,particle1List, particle2List, particle3List, particle4List, particle5List, particle6List, ...
        particle7List, particle8List, particle9List, particle10List, particle11List, particle12List, particle13List, particle14List, ...
        particle15List, particle16List, particle17List, particle18List ,particle19List, particle20List, concAgent, edgeChargeVector, ...
        particleTypeVector] = createBulkVector( g , porosity , concAgent , NZd , amountOfParticles , chargeAtEdges, edgeChargeVector );
    fileName    = ['InitialConfig/initialConfig','.', num2str(i),'.mat']; 
    save(fileName,'bulkVector','bulkTypeVector','chargeInfo','particle1List','particle2List','particle3List','particle4List','particle5List',...
        'particle6List','particle7List','particle8List','particle9List','particle10List','particle11List','particle12List','particle13List',...
        'particle14List','particle15List','particle16List','particle17List','particle18List','particle19List','particle20List','concAgent',...
        'edgeChargeVector','particleTypeVector')
    
end
