function [EW_dtens,pdist_matrix] = dtens_pdist_MainAggregatFolded(por)


%% Compiling the C++ Components, necessary to determine the particle size distribution
mex particleSizeDistribution.cpp;

%% Definition of parameters                                                 %% Can be changed
NX          = 25;                       % Number of horizontal Elements
NZd         = NX;                       % Number of vertical Elements
porosity    = por;                    % Porosity
numSmallPor = 0.05;                   % Approximate percentage of small pores
numBigPor   = 0.05;                   % Approximate percentage of big pores

singleBulkStencilLayers    = 3;          % Size of Stencil for Single Cell Bulk
doubleBulkStencilLayers    = 3;          % Size of Stencil for Double Bulk
trippleBulkStencilLayers   = 3;          % Size of Stencil for Tripple Bulk
quadrupleBulkStencilLayers = 3;          % Size of Stencil for Quadrupel Bulk
quintupleBulkStencilLayers = 3;          % Size of Stencil for Quintupel Bulk

tau         = 1;                        
epsilon     = 1;                        % Everything < epsilon is defined as zero
numInnerIt  = 10;                       % Number of Diffusion/Absorption Iterations
numOuterIt  = 20;                       % Number of Transport/Growth Iterations

p           = 1;                        % Order of Polynomials for LDG Disc.
ord         = 4;                        % Order of Gaussian Integration Rule
eta         = 1;                        % Penalty Parameter for LDG

startConBio = 0;                        % Initial Concentration of Biomass
maxConcBio  = 100;                      % Max. Concentration of Biomass -> Growth
minConcBio  = 10;                       % Min. Concentration of Biomass -> Shrinkage

f_uptakeAgent    = 1;                   % Growth Rate of Agent
f_decay          = -0.07;               % Degeneration Rate of Agent
f_uptakeAir      = -epsilon/2;          % Amount of consumed Air when Agent grows
f_bioMgrow       = 0; %10;                  % Growth Rate of Biomass
f_bioMdecr       = 0; %-0.1;                % Degeneration Rate of Biomass
f_bioUptAir      = 0; %-epsilon/3;          % Amount of Air consumed, when Biomass grows
f_uptakeBacteria = 0; %20;                  % Growth Rate of Bacteria
f_decayBac       = 0; %-epsilon/10;         % Degeneration Rate of Bacteria
f_oxyBac         = 0; %-epsilon/2;          % Amount of Air consumed, when Bacteria grows

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
f           = zero;
bacF        = zero;
bacG_N      = zero;
K11         = one;
K12         = zero;
K21         = zero;
K22         = one;
g_N         = zero;

%% Creating domain for Simulation
g = createDomainFolded(width, NX, NZd, NZu, intBound, upBound);
markEbdr        = g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8;
% concentration of sticky agent at t = 0
concAgent       = 10*ones( g.numE , 1 );    % Initial Concentration of Agent
NX = g.NX;

fileID = fopen( 'Move_bulk_log_file.odt' , 'w' );
%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed

% test

% bulkVector = singleCells( NX, 5 );
% bulkVector = doubleSingleCells( NX, 1 );
% bulkVector = chessgrid( NX );
% bulkVector = grid( NX, 5 );
% bulkVector = randi([0,10], g.numTsub, 1);
% bulkVector( find( bulkVector( : ) > 1 ) ) = 0;
% bioMvector = zeros(size(bulkVector));

helperPor = porosity;
while true
    bulkVector = createBulkVector(g, floor(NX*NX*helperPor*numBigPor), floor(NX*NX*helperPor*numSmallPor));
    if sum(~bulkVector) < NX^2 * (porosity - 0.01)
        helperPor = helperPor + 0.01;
    elseif sum(~bulkVector) > NX^2 * (porosity + 0.01)
        helperPor = helperPor - 0.01;
    else
        break
    end
end
bioMvector = zeros(size(bulkVector));


airVector = ~bulkVector;
bioMvec0airVec = bioMvector(airVector);
% K_hom = HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector)
% HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)))
% eigs(K)
numSolid = sum(bulkVector);

%% Defining Initial Concentrations of Biomass
bioConcVector   = startConBio * bioMvector;
%biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
%visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', 0);
bioMvector      = bioMvector > 0;

%% Computation of Local Matrices for LDG Scheme
computeBasesOnQuad(p, ord);

[hatMc, hatMx]          = computeHatM(p, ord);
[hatGc, hatGx, hatGy]   = computeHatG(p, ord);
[hatHc, hatHx, hatHy]   = computeHatH(p, ord);
hatRdiag                = computeHatRdiag(p, ord);
hatRoffdiag             = computeHatRoffdiag(p, ord);
hatSdiag                = computeHatSdiag(p, ord);
hatSoffdiag             = computeHatSoffdiag(p, ord);

%% Calculation of outer Loop for Movement of Bulk and Biomass Growth/Shrinkage
EW_dtens = zeros(2,numOuterIt); %List of EV of dtens
pdist_matrix    = zeros(30,2*numOuterIt); %matrix of the particleDistributions




for k = 1 : numOuterIt
    fprintf( fileID , 'step: %d \n' , k );
    
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
    globM       = assembleGlobMsub(g, hatMc, hatMx, airVectorSub, numAir);
    globH       = assembleGlobHsub(g, hatHc, hatHx, hatHy, airVectorSub, numAir);
    globQ       = assembleGlobQsub(g, markE0TairySub, hatSdiag, hatSoffdiag, airVectorSub, numAir);
    globQN      = assembleGlobQNsub(g, markE0TneumSub, hatSdiag, airVectorSub, numAir);
    globS       = assembleGlobSsub(g, markE0TairySub, hatSdiag, hatSoffdiag, eta, airVectorSub, numAir);

    if k == 1
        uDG         = projectAlg2DGsub(g, uZero, p, ord, hatMc, airVectorSub);
        q1DG        = projectAlg2DGsub(g, q1Zero, p, ord, hatMc, airVectorSub);
        q2DG        = projectAlg2DGsub(g, q2Zero, p, ord, hatMc, airVectorSub);
        uDGbac      = projectAlg2DGsub(g, bacZero, p, ord, hatMc, airVectorSub);
        q1DGbac     = projectAlg2DGsub(g, q1ZeroBac, p, ord, hatMc, airVectorSub);
        q2DGbac     = projectAlg2DGsub(g, q2ZeroBac, p, ord, hatMc, airVectorSub);
    else
        uDG         = uDG(airVectorSub,:);
        q1DG        = q1DG(airVectorSub,:);
        q2DG        = q2DG(airVectorSub,:);
        uDGbac      = uDGbac(airVectorSub,:);
        q1DGbac     = q1DGbac(airVectorSub,:);
        q2DGbac     = q2DGbac(airVectorSub,:);
    end  % if - else
    
    fDG         = projectAlg2DGsub(g, f, p, ord, hatMc, airVectorSub);
    fDGbac      = projectAlg2DGsub(g, bacF, p, ord, hatMc, airVectorSub);
    K11DG       = projectAlg2DGsub(g, K11, p, ord, hatMc, airVectorSub);
    K12DG       = projectAlg2DGsub(g, K12, p, ord, hatMc, airVectorSub);
    K21DG       = projectAlg2DGsub(g, K21, p, ord, hatMc, airVectorSub);
    K22DG       = projectAlg2DGsub(g, K22, p, ord, hatMc, airVectorSub);

    NT = numAir;
    N = (p+1)^2;

    globL       = globM * reshape(fDG', NT*N, 1);
    globLbac    = globM * reshape(fDGbac', NT*N, 1);
    sysU        = reshape(uDG', NT*N, 1);
    sysUbac     = reshape(uDGbac', NT*N, 1);
    sysQ1       = reshape(q1DG', NT*N, 1);
    sysQ1bac    = reshape(q1DGbac', NT*N, 1);
    sysQ2       = reshape(q2DG', NT*N, 1);
    sysQ2bac    = reshape(q2DGbac', NT*N, 1);

    globG       = assembleGlobGsub(g, hatGc, hatGx, hatGy, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
    globR       = assembleGlobRsub(g, markE0TairySub, hatRdiag, hatRoffdiag, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
    globKN      = assembleGlobKNsub(g, p, ord, markE0TneumSub, g_N, airVectorSub, numAir);
    globKNbac   = assembleGlobKNsub(g, p, ord, markE0TneumSub, bacG_N, airVectorSub, numAir);

    sysW = [    sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)                   ;
                sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)                   ;
                sparse(NT*N,NT*N)       ,   sparse(NT*N,NT*N)       ,   globM                               ];
    sysA = [    globM                   ,   sparse(NT*N,NT*N)       ,   globH{1} + globQ{1} + globQN{1}     ;
                sparse(NT*N,NT*N)       ,   globM                   ,   globH{2} + globQ{2} + globQN{2}     ;
                globG{1} + globR{1}     ,   globG{2} + globR{2}     ,   globS                               ];

    sysV = [    zeros(size(globM,1),1)  ;   zeros(size(globM,1),1)  ;   globKN + globL                      ];
    sysVbac = [ zeros(size(globM,1),1)  ;   zeros(size(globM,1),1)  ;   globKNbac + globLbac                ];
    sysX = [    sysQ1                   ;   sysQ2                   ;   sysU                                ];
    sysXbac = [  sysQ1bac                ;   sysQ2bac                ;   sysUbac                             ];
    
    if k == 1
        UDG    = zeros(g.numTsub, N);
        UDG(airVectorSub,:) = reshape( sysX( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
        UDGbac = zeros(g.numTsub, N);
        UDGbac(airVectorSub,:) = reshape( sysXbac( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
        uLagr       = projectDG2LagrangeSub( UDG );
        %bacLagr     = projectDG2LagrangeSub( UDGbac );
        visualizeDataSub(g, uLagr, 'u', 'solu', 0);
        %visualizeDataSub(g, bacLagr, 'bac', 'solBac', 0);
    end  % if

    %% Inner Loop for Diffusion and Growth/Shrinkage, ... (without Movement)
    
    for i = 1 : numInnerIt
        
        %% Calculations
        if size(uDG,1) == g.numTsub
            uDG = uDG(airVectorSub,:);
        end
        sysU        = reshape(uDG', NT*N, 1);
        sysX = [    sysQ1                   ;   sysQ2                   ;   sysU                                ];
        
        sysX = ( sysW + tau * sysA ) \ ( sysW * sysX + tau * sysV );
        sysXbac = ( sysW + tau * sysA ) \ ( sysW * sysXbac + tau * sysVbac );
        uDG         = zeros(g.numTsub, N);
        uDG(airVectorSub,:) = reshape( sysX( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
        uDGbac      = zeros(g.numTsub, N);
        uDGbac(airVectorSub,:) = reshape( sysXbac( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
%        uDG(~airVectorSub,1) = -1;
        
        concAir = computeConcAir(g, uDG, ord);
        concBac = computeConcAir(g, uDGbac, ord);

        %% Changes in Concentration of Agent       
        airyEdges = concAir(g.T0E(:,1)) > epsilon | concAir(g.T0E(:,2)) > epsilon;
        airyEdges = airyEdges & ( bulkVector(g.T0E(:,1)) | bulkVector(g.T0E(:,2)) | bioMvector(g.T0E(:,1)) | bioMvector(g.T0E(:,2)) );
        airyFactor = (concAir(g.T0E(:,1)) > epsilon) + (concAir(g.T0E(:,2)) > epsilon);
        airyTrapezoids = zeros(g.numT,1);
        airyTrapezoids(g.T0E(airyEdges, 1)) = 1;
        airyTrapezoids(g.T0E(airyEdges, 2)) = 1;
        airyTrapezoids = airyTrapezoids .* (concAir > epsilon);
        airyTrapezoids = airyTrapezoids > 0;
        concAgent(airyEdges) = concAgent(airyEdges) + tau * (f_uptakeAgent * airyFactor(airyEdges) + f_decay * concAgent(airyEdges));
        concAgent(~airyEdges) = concAgent(~airyEdges) + tau * f_decay * concAgent(~airyEdges);
        concAgent = concAgent .* (concAgent > epsilon);
        uDG(airyTrapezoids, 1) = uDG(airyTrapezoids, 1) + tau * f_uptakeAir;
        %% Changes in Concentration of Biomass
%         helperling = bioMvector .* (concBac > 0);
%         helperling = find(helperling);
%         bioConcVector(helperling) = bioConcVector(helperling) + concBac(helperling);
%         bioMvector = bioMvector | (concBac > minConcBio);
        bioConcVector(concBac > minConcBio) = concBac(concBac > minConcBio);
        bioMvector = bioConcVector > minConcBio;
        uDGbac((concBac > minConcBio),:) = 0;
        concBac(concBac > minConcBio) = 0;
        helperling = bioMvector .* (concBac > 0);
        helperling = find(helperling);
        helperlein = max(min(0.8, 1 - bioConcVector(helperling) / maxConcBio), 0.2) .* concBac(helperling);
        bioConcVector(helperling) = bioConcVector(helperling) + helperlein;
        uDGbac(helperling,1) = uDGbac(helperling,1) - helperlein;
        concBac(helperling) = concBac(helperling) - helperlein;
        bioMvector = bioConcVector > minConcBio;
%         helperlein = max(min(0.8, 1 - bioConcVector(bioMvector) / maxConcBio), 0.2) .* uDGbac(bioMvector);
%         helperlein2 = uDGbac(bioMvector,1) - helperlein;
%         uDGbac(bioMvector,1) = helperlein2;
%         bioConcVector(bioMvector) = bioConcVector(bioMvector) + helperlein;
%         bioConcVector(bioMvector) = max(bioConcVector(bioMvector), minConcBio);
        bioMair = bioMvector .* (concAir > epsilon);
        bioMair = bioMair > 0;
        bioMnoAir = bioMvector .* (concAir <= epsilon);
        bioMnoAir = bioMnoAir > 0;
        bioConcVector(bioMair) = bioConcVector(bioMair) + tau * f_bioMgrow;
%         bioConcVector(bioMnoAir) = bioConcVector(bioMnoAir) + tau * f_bioMdecr * bioConcVector(bioMnoAir);
% Implicit:
        bioConcVector(bioMnoAir) = 1/(1 - tau * f_bioMdecr) * bioConcVector(bioMnoAir);
        bioMvector = bioConcVector > minConcBio;
        bioConcVector = bioConcVector .* bioMvector;
        uDG(bioMair, 1) = uDG(bioMair, 1) + tau * f_bioUptAir;
        %% Changes in Concentration of Bacteria
        concAir = computeConcAir(g, uDG, ord);
        concAir2 = concAir > epsilon;
        concBac = computeConcAir(g, uDGbac, ord);
        uDGbac(concAir2,1) = uDGbac(concAir2,1) + tau * (f_uptakeBacteria + f_decayBac * uDGbac(concAir2,1));
        uDGbac(~concAir2,1) = uDGbac(~concAir2,1) + tau * f_decayBac * uDGbac(~concAir2,1);
        uDGbac(:,1) = uDGbac(:,1) .* (uDGbac(:,1) > 0);
        uDG(concAir2,1) = uDG(concAir2, 1) + tau * f_oxyBac;
        
%        concAir = computeConcAir(g, uDG, ord);
%        concAir = concAir < epsilon;
%        uDG(concAir,:) = 0;
    end  % for i

    %% Calculations
    q1DG        = zeros(g.numTsub, N);
    q1DG(airVectorSub,:) = reshape( sysX( 1 : NT * N ), N, NT )';
    q2DG        = zeros(g.numTsub, N);
    q2DG(airVectorSub,:) = reshape( sysX( NT * N + 1 : 2 * NT * N ), N, NT )';
    q1DGbac        = zeros(g.numTsub, N);
    q1DGbac(airVectorSub,:) = reshape( sysXbac( 1 : NT * N ), N, NT )';
    q2DGbac        = zeros(g.numTsub, N);
    q2DGbac(airVectorSub,:) = reshape( sysXbac( NT * N + 1 : 2 * NT * N ), N, NT )';
    
    %% Doing the Movement of Biomass
    airVector = ~bulkVector;
    bioOverHead = bioConcVector - maxConcBio;
    OverHead    = bioOverHead > 0;
    bioOverHead = bioOverHead .* OverHead;
%     bioConcVector(OverHead) = maxConcBio;
    
    numTrial = 0;
    
    while(numTrial < 10 && sum(OverHead) > 0)
        
    
    candidates = find(OverHead);
    %helper1 = [ candidates - NX , candidates - 1 , candidates , candidates + 1 , candidates + NX ];
    helper1 = [ NX * mod( floor( ( candidates - NX - 1 ) / NX ) , NZd ) + mod( candidates - NX - 1 , NX ) + 1 , ...
        NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 , ...
        candidates , ...
        NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 , ...
        NX * mod( floor( ( candidates + NX - 1 ) / NX ) , NZd ) + mod( candidates + NX - 1 , NX ) + 1 ];
    
    helper2 = (helper1 > 0) .* (helper1 < g.numT + 1);
    helper2 = helper2 > 0;
    aims = zeros(length(candidates),5);
    aims(helper2) = helper1(helper2);
    aims = reshape(aims, length(candidates), 5)';
    values = zeros(size(aims));
    
    for i = 1 : size(aims,1)
        for j = 1 : size(aims,2)
            if aims(i,j) > 0
                aims(i,j) = aims(i,j) .* airVector(aims(i,j)) .* ~OverHead(aims(i,j));
            end
        end
    end
    
    [maximo, indMax] = max(values);
    aims2 = zeros(size(indMax'));
    for i = 1 : length(indMax)
%         aims2(i) = aims(indMax(i),i);
      indices = find( maximo(i) == values(:,i) );
      chooser = randi(length(indices));
      if length(indices) == 1 && indMax(i) ~= indices(chooser)
          error('IDIOT')
      end
      aims2(i) = aims(indices(chooser),i);
    end
    aims = aims2;
    [n, bin] = histc(aims, unique(aims));
    multiple = find(n > 1);
    index    = find(ismember(bin, multiple));
    
    for i = 1 : length(index)-1
        if aims(index(i)) ~= 0
%         for j = i+1 : length(index)
%             if aims(index(i)) == aims(index(j))
%                 candidates(index(j)) = 0;
%                 aims(index(j)) = 0;
%             end  % if equal
%         end  % for j
            aims(index(i+1:end)) = aims(index(i+1:end)) .* (aims(index(i+1:end)) ~= aims(index(i)));
        end  % if zero
    end  % for i
    
    candidates = candidates(aims ~= 0);
    aims = aims(aims ~= 0);
    transport = min(bioOverHead(candidates), maxConcBio - bioConcVector(aims));
    bioConcVector(aims) = bioConcVector(aims) + transport;
    bioOverHead(candidates) = bioConcVector(candidates) - transport;
    OverHead    = bioOverHead > 0;
    
    
%     candidates = find(OverHead);
%     helper1 = [  candidates - 2 * NX , candidates - NX - 1 , candidates - NX , candidates - NX + 1 , ...
%             candidates - 2 , candidates - 1 , candidates , candidates + 1 , candidates + 2 , ...
%             candidates + NX - 1 , candidates + NX , candidates + NX + 1 , candidates + 2 * NX ];
    helper1 = [ NX * mod( floor( ( candidates - 2*NX - 1 ) / NX ) , NZd ) + mod( candidates - 2*NX - 1 , NX ) + 1 , ...
        NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) - NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) - NX - 1 , NX ) + 1 , ...
        NX * mod( floor( ( candidates - NX - 1 ) / NX ) , NZd ) + mod( candidates - NX - 1 , NX ) + 1 , ...
        NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 ) - NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 ) - NX - 1 , NX ) + 1 , ...
        NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 3, NX ) + 1 , ...
        NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 , ...
        candidates , ...
        NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 1 , ...
        NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + 2 , ...
        NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 , NX ) + 1 , ...
        NX * mod( floor( ( candidates + NX - 1 ) / NX ) , NZd ) + mod( candidates + NX - 1 , NX ) + 1 , ...
        NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates - 2, NX ) + 1 ) + NX - 1 , NX ) + 1 , ...
        NX * mod( floor( ( candidates + 2*NX - 1 ) / NX ) , NZd ) + mod( candidates + 2*NX - 1 , NX ) + 1 ];

    helper2 = (helper1 > 0) .* (helper1 < g.numT + 1);
    helper2 = helper2 > 0;
    aims = zeros(length(candidates),13);
    aims(helper2) = helper1(helper2);
    aims = reshape(aims, length(candidates), 13)';
    values = zeros(size(aims));
    
    for i = 1 : size(aims,1)
        for j = 1 : size(aims,2)
            if aims(i,j) > 0
                aims(i,j) = aims(i,j) .* airVector(aims(i,j)) .* ~OverHead(aims(i,j));
            end
        end
    end
    
    [maximo, indMax] = max(values);
    aims2 = zeros(size(indMax'));
    for i = 1 : length(indMax)
%         aims2(i) = aims(indMax(i),i);
      indices = find( maximo(i) == values(:,i) );
      chooser = randi(length(indices));
      if length(indices) == 1 && indMax(i) ~= indices(chooser)
          error('IDIOT')
      end
      aims2(i) = aims(indices(chooser),i);
    end
    aims = aims2;
    [n, bin] = histc(aims, unique(aims));
    multiple = find(n > 1);
    index    = find(ismember(bin, multiple));
    
    for i = 1 : length(index)-1
        if aims(index(i)) ~= 0
%         for j = i+1 : length(index)
%             if aims(index(i)) == aims(index(j))
%                 candidates(index(j)) = 0;
%                 aims(index(j)) = 0;
%             end  % if equal
%         end  % for j
            aims(index(i+1:end)) = aims(index(i+1:end)) .* (aims(index(i+1:end)) ~= aims(index(i)));
        end  % if zero
    end  % for i
    
    candidates = candidates(aims ~= 0);
    aims = aims(aims ~= 0);
    transport = min(bioOverHead(candidates), maxConcBio - bioConcVector(aims));
    bioConcVector(aims) = bioConcVector(aims) + transport;
    bioOverHead(candidates) = bioConcVector(candidates) - transport;
    OverHead    = bioOverHead > 0;
    numTrial = numTrial + 1;
    
    end % while sum(OverHead > 0)
    
    runNumber = 0;
    anz = sum(OverHead);
    while (runNumber < 10 && anz ~= 0)
        for i = 1 : anz
            culprit = find(OverHead);
            culprit = culprit(1);
            helpling = ~bulkVector .* ~OverHead;
            helpling = find(helpling);
            helpling2 = abs(helpling - culprit);
            distance = mod(helpling2, NX) + floor(helpling2/NX);
            
            [minimo, indMax] = min(distance);
            indices = find( minimo == distance );
            chooser = randi(length(indices));
            if length(indices) == 1 && indMax ~= indices(chooser)
                error('IDIOT')
            end
            transport = min(bioOverHead(culprit), maxConcBio - bioConcVector(helpling(indices(chooser))));
            bioConcVector(helpling(indices(chooser))) = bioConcVector(helpling(indices(chooser))) + transport;
            bioConcVector(culprit) = bioConcVector(culprit) - transport;            
        end  % for i
        
        bioOverHead = bioConcVector - maxConcBio;
        OverHead    = bioOverHead > 0;
        anz = sum(OverHead);
        runNumber = runNumber + 1;
    end  % while
    
    bioConcVector(OverHead) = maxConcBio;
    bioMvector = bioConcVector > minConcBio;
    bioConcVector = bioConcVector .* bioMvector;
    
    helperling = bioMvector .* (concBac > 0);
    helperling = find(helperling);
    bioMvector(helperling) = bioMvector(helperling) + concBac(helperling);
    uDGbac(bioMvector,:) = 0;

%% Doing the movement of Bulk
    
% Single Bulk

[bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent] = moveTripleBulk( 1, singleBulkStencilLayers, g, bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, markEbdr, N, bioMvector, NZd , fileID );

assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

% Double Bulk

[bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent] = moveTripleBulk( 2, doubleBulkStencilLayers, g, bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, markEbdr, N, bioMvector, NZd , fileID );

assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

% Triples of Bulk
 
[bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent] = moveTripleBulk( 3, trippleBulkStencilLayers, g, bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, markEbdr, N, bioMvector, NZd , fileID );

assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

% Quadruples of Bulk

[bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent] = moveTripleBulk( 4, quadrupleBulkStencilLayers, g, bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, markEbdr, N, bioMvector, NZd , fileID );

assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

% Quintuples of Bulk

[bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent] = moveTripleBulk( 5, quintupleBulkStencilLayers, g, bulkVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, concAgent, markEbdr, N, bioMvector, NZd , fileID );

assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nicht')

%     porosity = sum(airVector) / g.numT;
%     porosityBio = (sum(airVector) - sum(bioMvector)) / g.numT;

    
%     airVector = ~bulkVector;
%     bioMvec0airVec = bioMvector(airVector);
%     K_hom = HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector)
    dtensor = HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)));

    uLagr       = projectDG2LagrangeSub( uDG );
    visualizeDataSub(g, uLagr, 'u', 'solu', k);
    %bacLagr     = projectDG2LagrangeSub( uDGbac );
    %visualizeDataSub(g, bacLagr, 'bac', 'solBac', k);
    %biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
    %visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', k);

    particleDistribution = particleSizeDistribution(bulkVector.');
    numberOfSolidParticles = sum(particleDistribution(:,1) .* particleDistribution(:,2));
    

    if sum(isnan(dtensor))~=[0 0]
        EW_dtens(:,k)=[0;0];
    else
        EW_dtens(:,k)=eig(dtensor);
    end
    
    pdist_matrix(1:size(particleDistribution,1),2*k-1:2*k)   = particleDistribution;
 
    if sum(bioMvector .* bulkVector ~= 0)
        error('NEEEEEEEEEEEIN!')
    end
    
end  % for k
fclose( fileID );
end