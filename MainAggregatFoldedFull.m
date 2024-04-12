function MainAggregatFolded(numRun)

%% Compiling the C++ Components, necessary to determine the particle size distribution
mex particleSizeDistribution.cpp;

%% Definition of parameters                %% Can be changed
NX          = 250;                          % Number of horizontal Elements
NZd         = NX;                          % Number of vertical Elements
porosity    = 0.9;                         % Porosity

bilder = 1;                                % 0: nur Bild von Anfangs- und Endzustand, 1: Bild zu jedem Schritt, 2: Bild zu ersten 10 Zuständen und Endzustand, 3: Bild nach allen 100 Schritten

% IlliteCharge = 0;                   % Charge of Illite particles  plättchen 
% GoethitCharge = 0;                  % Charge of Goethit particles stäbchen 
% QuartzCharge = 0;                    % Charge of Quartz particles 

% QuartzRadius = 5;                          % Possible values: 5, 10, 20, 50
% GoethitLength = 5;                         % Possible values: 3, 5, 15
% IlliteSize = 6;                            % Possible values: 6, 24, 96

% agent = 0;                                 % sticky agent value

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
                 
amountOfParticles = [0,0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; %Amount of particle 1, particle 2, ... (sum should be 1)

attraction_type = 3; % 1: old volume charges, 2: no charges, 3: edge Charges

particle1StencilLayers = 5;
particle2StencilLayers = 2;
particle3StencilLayers = 6;
particle4StencilLayers = 3;
particle5StencilLayers = 5;
particle6StencilLayers = 6;
particle7StencilLayers = 6;
particle8StencilLayers = 3;
particle9StencilLayers = 4;
particle10StencilLayers = 4;
particle11StencilLayers = 4;
particle12StencilLayers = 4;
particle13StencilLayers = 3;
particle14StencilLayers = 3;
particle15StencilLayers = 3;
particle16StencilLayers = 3;
particle17StencilLayers = 3;
particle18StencilLayers = 3;
particle19StencilLayers = 3;
particle20StencilLayers = 3;


% QuartzBulkStencilLayers = 5;               % Sprungweite Quartz
% GoethitBulkStencilLayers = 5;              % Sprungweite Goethit
% IlliteBulkStencilLayers = 5;               % Sprungweite Illite

% bigParticleStencilLayersConstant = 4*particle3StencilLayers;  %<= 4*particle3StencilLayers
bigParticleStencilLayersConstant = 20;  %<= 4*particle3StencilLayers

% Amount = 0.20;
% str = pwd;
% Amount = str2double(['0.',str(end-6:end-5)]);



% tau         = 1;                        
% epsilon     = 1e-6;                     % Everything < epsilon is defined as zero
% numInnerIt  = 10;                       % Number of Diffusion/Absorption Iterations
numOuterIt  = 100;                     % Number of Transport/Growth Iterations

p           = 1;                        % Order of Polynomials for LDG Disc.
ord         = 4;                        % Order of Gaussian Integration Rule
% eta         = 1;                        % Penalty Parameter for LDG

startConBio = 0;                        % Initial Concentration of Biomass
% maxConcBio  = 100;                      % Max. Concentration of Biomass -> Growth
% minConcBio  = 10;                       % Min. Concentration of Biomass -> Shrinkage
% 
% f_uptakeAgent    = 0;                   % Growth Rate of Agent
% f_decay          = 0;                   % Degeneration Rate of Agent
% f_uptakeAir      = 0;                   % Amount of consumed Air when Agent grows
% f_bioMgrow       = 0; %10;              % Growth Rate of Biomass
% f_bioMdecr       = 0; %-0.1;            % Degeneration Rate of Biomass
% f_bioUptAir      = 0; %-epsilon/3;      % Amount of Air consumed, when Biomass grows
% f_uptakeBacteria = 0; %20;              % Growth Rate of Bacteria
% f_decayBac       = 0; %-epsilon/10;     % Degeneration Rate of Bacteria
% f_oxyBac         = 0; %-epsilon/2;      % Amount of Air consumed, when Bacteria grows

%% Fixed Parameters (of general framework)
NZu         = 0;
width       = NX;
intBound    = NZd * ones(NX+1,1);
upBound     = NZd * ones(NX+1,1);

zero        = @(x,y) x-x + 0;
one         = @(x,y) x-x + 1;

% concentration of nutrient u at t = 0
 uZero       = one;

% concentration of bacteria at t = 0
 bacZero     = zero;

q1Zero      = zero;
q2Zero      = zero;
q1ZeroBac   = zero;
q2ZeroBac   = zero;

% right hand side
% f           = zero;
% bacF        = zero;
% bacG_N      = zero;
% K11         = one;
% K12         = zero;
% K21         = zero;
% K22         = one;
% g_N         = zero;

%% Creating domain for Simulation
g = createDomainFolded(width, NX, NZd, NZu, intBound, upBound);
markEbdr        = g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8;
% concentration of sticky agent at t = 0
concAgent       = 0*ones( g.numE , 1 );    % Initial Concentration of Agent
edgeChargeVector = 0*ones( g.numCE , 1);    % Initial charge Vector
NX = g.NX;

fileID = fopen( 'Move_bulk_log_file' , 'w' );
% fileID_1 = fopen( 'Verteilungen', 'w' );
%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed

if nargin == 1
%     fileName    = ['initialConfig','.', num2str(numRun),'.mat'];
    load(['InitialConfig/initialConfig','.', num2str(numRun),'.mat'],'bulkVector','bulkTypeVector','chargeInfo','particle1List','particle2List',...
        'particle3List','particle4List','particle5List','particle6List','particle7List','particle8List','particle9List','particle10List',...
        'particle11List','particle12List','particle13List','particle14List','particle15List','particle16List','particle17List','particle18List',...
        'particle19List','particle20List','concAgent','edgeChargeVector','particleTypeVector');

else
    [bulkVector,bulkTypeVector,chargeInfo,particle1List, particle2List, particle3List, particle4List, particle5List, particle6List, particle7List,...
        particle8List, particle9List, particle10List, particle11List, particle12List, particle13List, particle14List, particle15List, ...
        particle16List, particle17List, particle18List, particle19List, particle20List, concAgent, edgeChargeVector, particleTypeVector] = ...
        createBulkVector( g , porosity , concAgent , NZd , amountOfParticles , chargeAtEdges, edgeChargeVector );
end
   
sumAgent = sum(concAgent);
bioMvector = zeros(size(bulkVector));

% airVector = ~bulkVector;
% bioMvec0airVec = bioMvector(airVector);
% K_hom = HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector)
% HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)))
% eigs(K)

numSolid = sum(bulkVector);

%% Defining Initial Concentrations of Biomass
% bioConcVector   = startConBio * bioMvector;
%biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
%visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', 0);
bioMvector      = bioMvector > 0;

%% Computation of Local Matrices for LDG Scheme
computeBasesOnQuad(p, ord);

[hatMc, hatMx]          = computeHatM(p, ord);
% [hatGc, hatGx, hatGy]   = computeHatG(p, ord);
% [hatHc, hatHx, hatHy]   = computeHatH(p, ord);
% hatRdiag                = computeHatRdiag(p, ord);
% hatRoffdiag             = computeHatRoffdiag(p, ord);
% hatSdiag                = computeHatSdiag(p, ord);
% hatSoffdiag             = computeHatSoffdiag(p, ord);

%% Calculation of outer Loop for Movement of Bulk and Biomass Growth/Shrinkage

for k = 1 : numOuterIt
    assert(sumAgent == sum(concAgent))
    fprintf('Step %d of %d\n', k, numOuterIt)
%     fprintf( fileID , 'step: %d \n' , k );
    
    assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
    
    %% Creating Vectors characterizing the Type of Edges
    idE                  = zeros(g.numE, 1);
    idE(g.E0T(bulkVector == 1,1)) = idE(g.E0T(bulkVector == 1,1)) + 1;
    idE(g.E0T(bulkVector == 1,2)) = idE(g.E0T(bulkVector == 1,2)) + 1;
    idE(g.E0T(bulkVector == 1,3)) = idE(g.E0T(bulkVector == 1,3)) + 1;
    idE(g.E0T(bulkVector == 1,4)) = idE(g.E0T(bulkVector == 1,4)) + 1;
    idE = idE + markEbdr;

    idEneum         = idE == 1;
    idEairy         = idE == 0;
    idEbulk         = ~idEairy;
    idEbulk2        = idE == 2;
    markE0Tneum     = idEneum(g.E0T);
    markE0Tairy     = idEairy(g.E0T);
    airVector   = ~bulkVector;

    markE0TneumSub = markE0Tneum(1:g.numTsub,:);
    markE0TairySub = markE0Tairy(1:g.numTsub,:);
    airVectorSub = airVector(1:g.numTsub);
    markE0TneumSub = markE0TneumSub(airVectorSub,:);
    markE0TairySub = markE0TairySub(airVectorSub,:);
    numAir = norm(airVectorSub.*ones(size(airVectorSub)),1);
    
    %% Starting PDE-time
    tic;

    %% Assembling global Matrices (independent of u, q, concentrations)
%     globM       = assembleGlobMsub(g, hatMc, hatMx, airVectorSub, numAir);
%     globH       = assembleGlobHsub(g, hatHc, hatHx, hatHy, airVectorSub, numAir);
%     globQ       = assembleGlobQsub(g, markE0TairySub, hatSdiag, hatSoffdiag, airVectorSub, numAir);
%     globQN      = assembleGlobQNsub(g, markE0TneumSub, hatSdiag, airVectorSub, numAir);
%     globS       = assembleGlobSsub(g, markE0TairySub, hatSdiag, hatSoffdiag, eta, airVectorSub, numAir);

    if k == 1
        uDG         = projectAlg2DGsub(g, uZero, p, ord, hatMc, airVectorSub);
        q1DG        = projectAlg2DGsub(g, q1Zero, p, ord, hatMc, airVectorSub);
        q2DG        = projectAlg2DGsub(g, q2Zero, p, ord, hatMc, airVectorSub);
%         uDGbac      = projectAlg2DGsub(g, bacZero, p, ord, hatMc, airVectorSub);
%         q1DGbac     = projectAlg2DGsub(g, q1ZeroBac, p, ord, hatMc, airVectorSub);
%         q2DGbac     = projectAlg2DGsub(g, q2ZeroBac, p, ord, hatMc, airVectorSub);
    else
        uDG         = uDG(airVectorSub,:);
        q1DG        = q1DG(airVectorSub,:);
        q2DG        = q2DG(airVectorSub,:);
%         uDGbac      = uDGbac(airVectorSub,:);
%         q1DGbac     = q1DGbac(airVectorSub,:);
%         q2DGbac     = q2DGbac(airVectorSub,:);
    end  % if - else
%     
%     fDG         = projectAlg2DGsub(g, f, p, ord, hatMc, airVectorSub);
%     fDGbac      = projectAlg2DGsub(g, bacF, p, ord, hatMc, airVectorSub);
%     K11DG       = projectAlg2DGsub(g, K11, p, ord, hatMc, airVectorSub);
%     K12DG       = projectAlg2DGsub(g, K12, p, ord, hatMc, airVectorSub);
%     K21DG       = projectAlg2DGsub(g, K21, p, ord, hatMc, airVectorSub);
%     K22DG       = projectAlg2DGsub(g, K22, p, ord, hatMc, airVectorSub);

    NT = numAir;
    N = (p+1)^2;

%     globL       = globM * reshape(fDG', NT*N, 1);
%     globLbac    = globM * reshape(fDGbac', NT*N, 1);
    sysU        = reshape(uDG', NT*N, 1);
%     sysUbac     = reshape(uDGbac', NT*N, 1);
    sysQ1       = reshape(q1DG', NT*N, 1);
%     sysQ1bac    = reshape(q1DGbac', NT*N, 1);
    sysQ2       = reshape(q2DG', NT*N, 1);
%     sysQ2bac    = reshape(q2DGbac', NT*N, 1);
% 
%     globG       = assembleGlobGsub(g, hatGc, hatGx, hatGy, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
%     globR       = assembleGlobRsub(g, markE0TairySub, hatRdiag, hatRoffdiag, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
%     globKN      = assembleGlobKNsub(g, p, ord, markE0TneumSub, g_N, airVectorSub, numAir);
%     globKNbac   = assembleGlobKNsub(g, p, ord, markE0TneumSub, bacG_N, airVectorSub, numAir);
% 
%     sysW = [    sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)                   ;
%                 sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)                   ;
%                 sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)       ,   globM                               ];
%     sysA = [    globM                   ,   sparse(NT*N,NT*N)       ,   globH{1} + globQ{1} + globQN{1}     ;
%                 sparse(NT*N,NT*N)       ,   globM                   ,   globH{2} + globQ{2} + globQN{2}     ;
%                 globG{1} + globR{1}     ,   globG{2} + globR{2}     ,   globS                               ];
% 
%     sysV = [    zeros(size(globM,1),1)  ;   zeros(size(globM,1),1)  ;   globKN + globL                      ];
%     sysVbac = [ zeros(size(globM,1),1)  ;   zeros(size(globM,1),1)  ;   globKNbac + globLbac                ];
    sysX = [    sysQ1                   ;   sysQ2                   ;   sysU                                ];
%     sysXbac = [  sysQ1bac                ;   sysQ2bac                ;   sysUbac                             ];
    
    if k == 1
        UDG    = zeros(g.numTsub, N);
        UDG(airVectorSub,:) = reshape( sysX( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
%         UDG(bulkTypeVector == QuartzCharge,1) = -1;
%         UDG(bulkVector == 1,1) = -1;
        q1DG        = UDG;
        q2DG        = UDG;
        uDGbac      = UDG;
        q1DGbac     = UDG;
        q2DGbac     = UDG;
        UDG(bulkTypeVector == -1,1) = -1;
        UDG(bulkTypeVector == 1,1) = 0;
%         UDGbac = zeros(g.numTsub, N);
%         UDGbac(airVectorSub,:) = reshape( sysXbac( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
%         uLagr       = projectDG2LagrangeSub( UDG );
        %bacLagr     = projectDG2LagrangeSub( UDGbac );
%         visualizeDataSub(g, uLagr, 'u', 'solu', 0);
        visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', 0);
        uDG         = UDG;
        printInfo(0,bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List, particle6List, particle7List,...
            particle8List,particle9List,particle10List,particle11List,particle12List,particle13List,particle14List,particle15List,...
            particle16List,particle17List,particle18List,particle19List,particle20List);
    else
        UDG    = zeros(g.numTsub, N);
        UDG(airVectorSub,:) = reshape( sysX( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
%         UDG(bulkTypeVector == QuartzCharge,1) = -1;
%         UDG(bulkVector == 1,1) = -1;
        q1DG        = UDG;
        q2DG        = UDG;
        uDGbac      = UDG;
        q1DGbac     = UDG;
        q2DGbac     = UDG; 
        UDG(bulkTypeVector == -1,1) = -1;
        UDG(bulkTypeVector == 1,1) = 0;
%         UDGbac = zeros(g.numTsub, N);
%         UDGbac(airVectorSub,:) = reshape( sysXbac( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
        %bacLagr     = projectDG2LagrangeSub( UDGbac );
        uDG         = UDG;  
        %visualizeDataSub(g, bacLagr, 'bac', 'solBac', 0);
    end  % if

    %% Inner Loop for Diffusion and Growth/Shrinkage, ... (without Movement)
    
%     for i = 1 : numInnerIt
%         
%         %% Calculations
%         if size(uDG,1) == g.numTsub
%             uDG = uDG(airVectorSub,:);
%         end
%         sysU        = reshape(uDG', NT*N, 1);
%         sysX = [    sysQ1                   ;   sysQ2                   ;   sysU                                ];
% %         
% %         sysX = ( sysW + tau * sysA ) \ ( sysW * sysX + tau * sysV );
% %         sysXbac = ( sysW + tau * sysA ) \ ( sysW * sysXbac + tau * sysVbac );
%         uDG         = zeros(g.numTsub, N);
%         uDG(airVectorSub,:) = reshape( sysX( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
% % 
% %         uDGbac      = zeros(g.numTsub, N);
% %         uDGbac(airVectorSub,:) = reshape( sysXbac( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
% % %        uDG(~airVectorSub,1) = -1;
% %         
%         concAir = computeConcAir(g, uDG, ord);
%         concBac = computeConcAir(g, uDGbac, ord);
% 
%         %% Changes in Concentration of Agent       
%         airyEdges = concAir(g.T0E(:,1)) > epsilon | concAir(g.T0E(:,2)) > epsilon;
%         airyEdges = airyEdges & ( bulkVector(g.T0E(:,1)) | bulkVector(g.T0E(:,2)) | bioMvector(g.T0E(:,1)) | bioMvector(g.T0E(:,2)) );
%         airyFactor = (concAir(g.T0E(:,1)) > epsilon) + (concAir(g.T0E(:,2)) > epsilon);
%         airyTrapezoids = zeros(g.numT,1);
%         airyTrapezoids(g.T0E(airyEdges, 1)) = 1;
%         airyTrapezoids(g.T0E(airyEdges, 2)) = 1;
%         airyTrapezoids = airyTrapezoids .* (concAir > epsilon);
%         airyTrapezoids = airyTrapezoids > 0;
%         assert(sumAgent == sum(concAgent))
%         concAgent(airyEdges) = concAgent(airyEdges) + tau * (f_uptakeAgent * airyFactor(airyEdges) + f_decay * concAgent(airyEdges));
%         assert(sumAgent == sum(concAgent))
%         concAgent(~airyEdges) = concAgent(~airyEdges) + tau * f_decay * concAgent(~airyEdges);
%         assert(sumAgent == sum(concAgent))
% %         assert(min(concAgent == 0))
%         concAgent = concAgent .* (concAgent > epsilon);
%         assert(sumAgent == sum(concAgent))
%         uDG(airyTrapezoids, 1) = uDG(airyTrapezoids, 1) + tau * f_uptakeAir;
%         %% Changes in Concentration of Biomass
% %         helperling = bioMvector .* (concBac > 0);
% %         helperling = find(helperling);
% %         bioConcVector(helperling) = bioConcVector(helperling) + concBac(helperling);
% %         bioMvector = bioMvector | (concBac > minConcBio);
%         bioConcVector(concBac > minConcBio) = concBac(concBac > minConcBio);
%         bioMvector = bioConcVector > minConcBio;
%         uDGbac((concBac > minConcBio),:) = 0;
%         concBac(concBac > minConcBio) = 0;
%         helperling = bioMvector .* (concBac > 0);
%         helperling = find(helperling);
%         helperlein = max(min(0.8, 1 - bioConcVector(helperling) / maxConcBio), 0.2) .* concBac(helperling);
%         bioConcVector(helperling) = bioConcVector(helperling) + helperlein;
%         uDGbac(helperling,1) = uDGbac(helperling,1) - helperlein;
%         concBac(helperling) = concBac(helperling) - helperlein;
%         bioMvector = bioConcVector > minConcBio;
% %         helperlein = max(min(0.8, 1 - bioConcVector(bioMvector) / maxConcBio), 0.2) .* uDGbac(bioMvector);
% %         helperlein2 = uDGbac(bioMvector,1) - helperlein;
% %         uDGbac(bioMvector,1) = helperlein2;
% %         bioConcVector(bioMvector) = bioConcVector(bioMvector) + helperlein;
% %         bioConcVector(bioMvector) = max(bioConcVector(bioMvector), minConcBio);
%         bioMair = bioMvector .* (concAir > epsilon);
%         bioMair = bioMair > 0;
%         bioMnoAir = bioMvector .* (concAir <= epsilon);
%         bioMnoAir = bioMnoAir > 0;
%         bioConcVector(bioMair) = bioConcVector(bioMair) + tau * f_bioMgrow;
% %         bioConcVector(bioMnoAir) = bioConcVector(bioMnoAir) + tau * f_bioMdecr * bioConcVector(bioMnoAir);
% % Implicit:
%         bioConcVector(bioMnoAir) = 1/(1 - tau * f_bioMdecr) * bioConcVector(bioMnoAir);
%         bioMvector = bioConcVector > minConcBio;
%         bioConcVector = bioConcVector .* bioMvector;
%         uDG(bioMair, 1) = uDG(bioMair, 1) + tau * f_bioUptAir;
%         %% Changes in Concentration of Bacteria
%         concAir = computeConcAir(g, uDG, ord);
%         concAir2 = concAir > epsilon;
%         concBac = computeConcAir(g, uDGbac, ord);
%         uDGbac(concAir2,1) = uDGbac(concAir2,1) + tau * (f_uptakeBacteria + f_decayBac * uDGbac(concAir2,1));
%         uDGbac(~concAir2,1) = uDGbac(~concAir2,1) + tau * f_decayBac * uDGbac(~concAir2,1);
%         uDGbac(:,1) = uDGbac(:,1) .* (uDGbac(:,1) > 0);
%         uDG(concAir2,1) = uDG(concAir2, 1) + tau * f_oxyBac;
% 
% %        concAir = computeConcAir(g, uDG, ord);
% %        concAir = concAir < epsilon;
% %        uDG(concAir,:) = 0;
%     end  % for i

    %% Calculations
%     q1DG        = zeros(g.numTsub, N);
%     q1DG(airVectorSub,:) = reshape( sysX( 1 : NT * N ), N, NT )';
%     q2DG        = zeros(g.numTsub, N);
%     q2DG(airVectorSub,:) = reshape( sysX( NT * N + 1 : 2 * NT * N ), N, NT )';
%     q1DGbac        = zeros(g.numTsub, N);
%     q1DGbac(airVectorSub,:) = reshape( sysXbac( 1 : NT * N ), N, NT )';
%     q2DGbac        = zeros(g.numTsub, N);
%     q2DGbac(airVectorSub,:) = reshape( sysXbac( NT * N + 1 : 2 * NT * N ), N, NT )';
    
    %% Doing the Movement of Biomass
%     airVector = ~bulkVector;
%     bioOverHead = bioConcVector - maxConcBio;
%     OverHead    = bioOverHead > 0;
%     bioOverHead = bioOverHead .* OverHead;
% %     bioConcVector(OverHead) = maxConcBio;
%     
%     numTrial = 0;
%     
%     while(numTrial < 10 && sum(OverHead) > 0)
%         
%     
%     candidates = find(OverHead);
%     %helper1 = [ candidates - NX , candidates - 1 , candidates , candidates + 1 , candidates + NX ];
%     helper1 = [ NX * mod( floor( ( candidates - NX - 1 ) / NX ) , NZd ) + mod( candidates - NX - 1 , NX ) + 1 , ...
%         NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 , ...
%         candidates , ...
%         NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 , ...
%         NX * mod( floor( ( candidates + NX - 1 ) / NX ) , NZd ) + mod( candidates + NX - 1 , NX ) + 1 ];
%     
%     helper2 = (helper1 > 0) .* (helper1 < g.numT + 1);
%     helper2 = helper2 > 0;
%     aims = zeros(length(candidates),5);
%     aims(helper2) = helper1(helper2);
%     aims = reshape(aims, length(candidates), 5)';
%     values = zeros(size(aims));
%     
%     for i = 1 : size(aims,1)
%         for j = 1 : size(aims,2)
%             if aims(i,j) > 0
%                 aims(i,j) = aims(i,j) .* airVector(aims(i,j)) .* ~OverHead(aims(i,j));Ftoc
%             end
%         end
%     end
%     
%     [maximo, indMax] = max(values);
%     aims2 = zeros(size(indMax'));
%     for i = 1 : length(indMax)
% %         aims2(i) = aims(indMax(i),i);
%       indices = find( maximo(i) == values(:,i) );
%       chooser = randi(length(indices));
%       if length(indices) == 1 && indMax(i) ~= indices(chooser)
%           error('IDIOT')
%       end
%       aims2(i) = aims(indices(chooser),i);
%     end
%     aims = aims2;
%     [n, bin] = histc(aims, unique(aims));
%     multiple = find(n > 1);
%     index    = find(ismember(bin, multiple));
%     
%     for i = 1 : length(index)-1
%         if aims(index(i)) ~= 0
% %         for j = i+1 : length(index)
% %             if aims(index(i)) == aims(index(j))
% %                 candidates(index(j)) = 0;
% %                 aims(index(j)) = 0;
% %             end  % if equal
% %         end  % for j
%             aims(index(i+1:end)) = aims(index(i+1:end)) .* (aims(index(i+1:end)) ~= aims(index(i)));
%         end  % if zero
%     end  % for i
%     
%     candidates = candidates(aims ~= 0);
%     aims = aims(aims ~= 0);
%     transport = min(bioOverHead(candidates), maxConcBio - bioConcVector(aims));
%     bioConcVector(aims) = bioConcVector(aims) + transport;
%     bioOverHead(candidates) = bioConcVector(candidates) - transport;
%     OverHead    = bioOverHead > 0;
%     
%     
% %     candidates = find(OverHead);
% %     helper1 = [  candidates - 2 * NX , candidates - NX - 1 , candidates - NX , candidates - NX + 1 , ...
% %             candidates - 2 , candidates - 1 , candidates , candidates + 1 , candidates + 2 , ...
% %             candidates + NX - 1 , candidates + NX , candidates + NX + 1 , candidates + 2 * NX ];
%     helper1 = [ NX * mod( floor( ( candidates - 2*NX - 1 ) / NX ) , NZd ) + mod( candidates - 2*NX - 1 , NX ) + 1 , ...
%         NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) - NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) - NX - 1 , NX ) + 1 , ...
%         NX * mod( floor( ( candidates - NX - 1 ) / NX ) , NZd ) + mod( candidates - NX - 1 , NX ) + 1 , ...
%         NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 ) - NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 ) - NX - 1 , NX ) + 1 , ...
%         NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 3, NX ) + 1 , ...
%         NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 , ...
%         candidates , ...
%         NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 , ...
%         NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 2 , ...
%         NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 , NX ) + 1 , ...
%         NX * mod( floor( ( candidates + NX - 1 ) / NX ) , NZd ) + mod( candidates + NX - 1 , NX ) + 1 , ...
%         NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 , NX ) + 1 , ...
%         NX * mod( floor( ( candidates + 2*NX - 1 ) / NX ) , NZd ) + mod( candidates + 2*NX - 1 , NX ) + 1 ];
% 
%     helper2 = (helper1 > 0) .* (helper1 < g.numT + 1);
%     helper2 = helper2 > 0;
%     aims = zeros(length(candidates),13);
%     aims(helper2) = helper1(helper2);
%     aims = reshape(aims, length(candidates), 13)';
%     values = zeros(size(aims));
%     
%     for i = 1 : size(aims,1)
%         for j = 1 : size(aims,2)
%             if aims(i,j) > 0
%                 aims(i,j) = aims(i,j) .* airVector(aims(i,j)) .* ~OverHead(aims(i,j));
%             end
%         end
%     end
%     
%     [maximo, indMax] = max(values);
%     aims2 = zeros(size(indMax'));
%     for i = 1 : length(indMax)
% %         aims2(i) = aims(indMax(i),i);
%       indices = find( maximo(i) == values(:,i) );
%       chooser = randi(length(indices));
%       if length(indices) == 1 && indMax(i) ~= indices(chooser)
%           error('IDIOT')
%       end
%       aims2(i) = aims(indices(chooser),i);
%     end
%     aims = aims2;
%     [n, bin] = histc(aims, unique(aims));
%     multiple = find(n > 1);
%     index    = find(ismember(bin, multiple));
%     
%     for i = 1 : length(index)-1
%         if aims(index(i)) ~= 0
% %         for j = i+1 : length(index)
% %             if aims(index(i)) == aims(index(j))
% %                 candidates(index(j)) = 0;
% %                 aims(index(j)) = 0;
% %             end  % if equal
% %         end  % for j
%             aims(index(i+1:end)) = aims(index(i+1:end)) .* (aims(index(i+1:end)) ~= aims(index(i)));
%         end  % if zero
%     end  % for i
%     
%     candidates = candidates(aims ~= 0);
%     aims = aims(aims ~= 0);
%     transport = min(bioOverHead(candidates), maxConcBio - bioConcVector(aims));
%     bioConcVector(aims) = bioConcVector(aims) + transport;
%     bioOverHead(candidates) = bioConcVector(candidates) - transport;
%     OverHead    = bioOverHead > 0;
%     numTrial = numTrial + 1;
%     
%     end % while sum(OverHead > 0)
%     
%     runNumber = 0;
%     anz = sum(OverHead);
%     while (runNumber < 10 && anz ~= 0)
%         for i = 1 : anz
%             culprit = find(OverHead);
%             culprit = culprit(1);
%             helpling = ~bulkVector .* ~OverHead;
%             helpling = find(helpling);
%             helpling2 = abs(helpling - culprit);
%             distance = mod(helpling2, NX) + floor(helpling2/NX);
%             
%             [minimo, indMax] = min(distance);
%             indices = find( minimo == distance );
%             chooser = randi(length(indices));
%             if length(indices) == 1 && indMax ~= indices(chooser)
%                 error('IDIOT')
%             end
%             transport = min(bioOverHead(culprit), maxConcBio - bioConcVector(helpling(indices(chooser))));
%             bioConcVector(helpling(indices(chooser))) = bioConcVector(helpling(indices(chooser))) + transport;
%             bioConcVector(culprit) = bioConcVector(culprit) - transport;            
%         end  % for i
%         
%         bioOverHead = bioConcVector - maxConcBio;
%         OverHead    = bioOverHead > 0;
%         anz = sum(OverHead);
%         runNumber = runNumber + 1;
%     end  % while
%     
%     bioConcVector(OverHead) = maxConcBio;
%     bioMvector = bioConcVector > minConcBio;
%     bioConcVector = bioConcVector .* bioMvector;
%     
%     helperling = bioMvector .* (concBac > 0);
%     helperling = find(helperling);
%     bioMvector(helperling) = bioMvector(helperling) + concBac(helperling);
%     uDGbac(bioMvector,:) = 0;
% 
     fprintf('Time of unnecessary stuff: %d ', toc)

%% Doing the movement of Bulk
assert(sumAgent == sum(concAgent),'sumAgent')

assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

% % Quartz 
% if QuartzRadius == 5
%     load('circle_5.mat','circle');
% elseif QuartzRadius == 10
%     load('circle_10.mat','circle');
% elseif QuartzRadius == 20
%     load('circle_20.mat','circle');
% elseif QuartzRadius == 50
%     load('circle_50.mat','circle');
% else
%     assert('wrong Quartzradius')
% end
% %% Quartz, typeflag = 1
% assert(sumAgent==sum(concAgent),'sumAgent')
% [bulkVector,bulkTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, QuartzList,~] = moveTripleBulk_ladung( length(circle), QuartzBulkStencilLayers, g, bulkVector, bulkTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, agent, markEbdr, N, bioMvector, NZd , fileID ,QuartzList,sumAgent,1,0);
% assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

tic

%% 
% particle_index = 1 for 1x17 (Goethit), 2 for 2x34 (Goethit), 3 for 2x6
% (Illite), 4 for 2x30  (Illite), 5 for 8x100 (Illite)

%% particle 1 (1x17, Goethit)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle1List,~] = moveTripleBulk_ladung( 17, particle1StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle1List,sumAgent,3,0, attraction_type,1);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 2 (2x34, Goethit)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle2List,~] = moveTripleBulk_ladung( 68, particle2StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle2List,sumAgent,3,0, attraction_type,2);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 3 (2x6, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle3List,~] = moveTripleBulk_ladung( 12, particle3StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle3List,sumAgent,3,0, attraction_type,3);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 4 (2x30, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle4List,~] = moveTripleBulk_ladung( 60, particle4StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle4List,sumAgent,3,0, attraction_type,4);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 5 (8x100, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle5List,~] = moveTripleBulk_ladung( 800, particle5StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle5List,sumAgent,3,0, attraction_type,5);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 6 (1x12, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle6List,~] = moveTripleBulk_ladung( 12, particle6StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle6List,sumAgent,3,0, attraction_type,6);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 7 (3x4, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle7List,~] = moveTripleBulk_ladung( 12, particle7StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle7List,sumAgent,3,0, attraction_type,7);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 8 (1x36, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle8List,~] = moveTripleBulk_ladung( 36, particle8StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle8List,sumAgent,3,0, attraction_type,8);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 9 (1x30, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle9List,~] = moveTripleBulk_ladung( 30, particle9StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle9List,sumAgent,3,0, attraction_type,9);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 10 (2x15, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle10List,~] = moveTripleBulk_ladung( 30, particle10StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle10List,sumAgent,3,0, attraction_type,10);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 11 (3x12, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle11List,~] = moveTripleBulk_ladung( 30, particle11StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle11List,sumAgent,3,0, attraction_type,11);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 12 (5x6, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle12List,~] = moveTripleBulk_ladung( 60, particle12StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle12List,sumAgent,3,0, attraction_type,12);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 13 (1x60, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle13List,~] = moveTripleBulk_ladung( 60, particle13StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle13List,sumAgent,3,0, attraction_type,13);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 14 (3x20, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle14List,~] = moveTripleBulk_ladung( 60, particle14StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle14List,sumAgent,3,0, attraction_type,14);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 15 (4x15, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle15List,~] = moveTripleBulk_ladung( 60, particle15StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle15List,sumAgent,3,0, attraction_type,15);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 16 (5x12, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle16List,~] = moveTripleBulk_ladung( 60, particle16StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle16List,sumAgent,3,0, attraction_type,16);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 17 (1x48, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle17List,~] = moveTripleBulk_ladung( 48, particle17StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle17List,sumAgent,3,0, attraction_type,17);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 18 (2x24, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle18List,~] = moveTripleBulk_ladung( 48, particle18StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle18List,sumAgent,3,0, attraction_type,18);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 19 (3x16, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle19List,~] = moveTripleBulk_ladung( 48, particle19StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle19List,sumAgent,3,0, attraction_type,19);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')
%% particle 20 (6x8, Illite)
assert(sumAgent==sum(concAgent),'sumAgent')
[bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, particle20List,~] = moveTripleBulk_ladung( 48, particle20StencilLayers, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particle20List,sumAgent,3,0, attraction_type,20);
assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')


fprintf('Time for jumping of single particles: %d ', toc)


% %% Goethit, typeflag = 2 (Stäbchen)
% assert(sumAgent==sum(concAgent),'sumAgent')
% [bulkVector,bulkTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, GoethitList,~] = moveTripleBulk_ladung( GoethitLength, GoethitBulkStencilLayers, g, bulkVector, bulkTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, agent, markEbdr, N, bioMvector, NZd , fileID ,GoethitList,sumAgent,2,0);
% assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

% %% Illite, typeflag = 3 (Plättchen)
% assert(sumAgent==sum(concAgent),'sumAgent')
% [bulkVector,bulkTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, IlliteList,~] = moveTripleBulk_ladung( IlliteSize, IlliteBulkStencilLayers, g, bulkVector, bulkTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, agent, markEbdr, N, bioMvector, NZd , fileID ,IlliteList,sumAgent,3,0);
% assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

%%
%% move big particles, typeflag = 4
tic
bigJumping = 1;
if bigJumping == 1
[particleList, particleContent] = particleInfoNew(bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List,particle6List,...
    particle7List,particle8List,particle9List,particle10List,particle11List,particle12List,particle13List,particle14List,particle15List,...
    particle16List,particle17List,particle18List,particle19List,particle20List);
for particle = 1 : length( particleList )
    particleSize = length( particleList{ particle } ); 
    if(size(particleContent{particle},1)<2)% kein Verbund
        continue
    end
    test_ind = particleList{particle}(1);
    
    bigParticleStencilLayers_individual = round(bigParticleStencilLayersConstant*1/(particleSize)^0.5);
    
%     if bigParticleStencilLayers_individual > QuartzBulkStencilLayers
%        bigParticleStencilLayers_individual = QuartzBulkStencilLayers;
%     end
    if bigParticleStencilLayers_individual == 0
        continue
    end
    assert(sumAgent==sum(concAgent),'sumAgent')
    [bulkVector,bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, changedList,~] = moveTripleBulk_ladung( particleSize, bigParticleStencilLayers_individual, g, bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, edgeChargeVector, 0, markEbdr, N, bioMvector, NZd , fileID ,particleList{ particle },sumAgent,4,0, attraction_type);
    assert(sumAgent==sum(concAgent),'sumAgent')
  
    % switch inds in QuartzLits, GoethitList, IlliteList
    if test_ind ~= changedList(1) % Listen müssen angepasst werden
        sten_ind = find(stencil(g.NX,g.NX,test_ind,bigParticleStencilLayers_individual)==changedList(1));
        for part = 1: size(particleContent{particle},1)  
            switch particleContent{particle}(part,1)
                case 1 
                    sten = stencil(g.NX,g.NX,particle1List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle1List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
                case 2 
                    sten = stencil(g.NX,g.NX,particle2List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle2List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
                case 3 
                    sten = stencil(g.NX,g.NX,particle3List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle3List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
                case 4 
                    sten = stencil(g.NX,g.NX,particle4List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle4List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';                   
                case 5
                    sten = stencil(g.NX,g.NX,particle5List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle5List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
                case 6
                    sten = stencil(g.NX,g.NX,particle6List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle6List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
               case 7
                    sten = stencil(g.NX,g.NX,particle7List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle7List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
               case 8
                    sten = stencil(g.NX,g.NX,particle8List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle8List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
               case 9
                    sten = stencil(g.NX,g.NX,particle9List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle9List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
               case 10
                    sten = stencil(g.NX,g.NX,particle10List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle10List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
              case 11
                    sten = stencil(g.NX,g.NX,particle11List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle11List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';
              case 12
                    sten = stencil(g.NX,g.NX,particle12List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle12List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 13
                    sten = stencil(g.NX,g.NX,particle13List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle13List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 14
                    sten = stencil(g.NX,g.NX,particle14List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle14List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 15
                    sten = stencil(g.NX,g.NX,particle15List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle15List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 16
                    sten = stencil(g.NX,g.NX,particle16List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle16List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 17
                    sten = stencil(g.NX,g.NX,particle17List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle17List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 18
                    sten = stencil(g.NX,g.NX,particle18List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle18List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 19
                    sten = stencil(g.NX,g.NX,particle19List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle19List(particleContent{particle}(part,2),:)=sten(:,sten_ind)'; 
              case 20
                    sten = stencil(g.NX,g.NX,particle20List( particleContent{particle}(part,2) ,:),bigParticleStencilLayers_individual);
                    particle20List(particleContent{particle}(part,2),:)=sten(:,sten_ind)';       
            end
        end
    end
    assert(sumAgent==sum(concAgent),'sumAgent_rotation_before')
end
end
fprintf('Time for jumping of aggregates: %d ', toc)
%%
%     porosity = sum(airVector) / g.numT;
%     porosityBio = (sum(airVector) - sum(bioMvector)) / g.numT;

    
%     airVector = ~bulkVector;
%     bioMvec0airVec = bioMvector(airVector);
%     K_hom = HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector)
%     dtensor=HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)))

%     uDG(bulkTypeVector == QuartzCharge,1) = -1;
%     uDG(bulkVector == 1,1) = -1;
    
    uDG(bulkTypeVector == -1,1) = -1;
    uDG(bulkTypeVector == 1,1) = 0;
    
    %alt
    if bilder == 0 && k == numOuterIt
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k);
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k);
    elseif bilder == 1
%     uLagr       = projectDG2LagrangeSub( uDG );
%     visualizeDataSub(g, uLagr, 'u', 'solu', k);  
    visualizeDataSub(g, particleTypeVector, 'particleType', 'solu', k); 
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
    


printInfo(k,bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List, particle6List, particle7List, particle8List, ...
    particle9List, particle10List, particle11List, particle12List, particle13List, particle14List, particle15List, particle16List, particle17List, ...
    particle18List, particle19List, particle20List)

%     numberOfSolidParticles = sum(particleDistribution(:,1) .* particleDistribution(:,2));

    if sum(bioMvector .* bulkVector ~= 0)
        error('NEEEEEEEEEEEIN!')
    end
%     particleDistribution

    if(mod(k,1000) == 0 || k == numOuterIt)
        fileName    = ['config','.', num2str(k),'.mat']; 
        save(fileName,'bulkVector','bulkTypeVector','chargeInfo','particle1List','particle2List','particle3List','particle4List','particle5List',...
            'particle6List','particle7List','particle8List','particle9List','particle10List','particle11List','particle12List','particle13List',...
            'particle14List','particle15List','particle16List','particle17List','particle18List','particle19List','particle20List','concAgent','edgeChargeVector')        
    end
    
end  % for k

fclose( fileID );
% fclose( fileID_1);
end
