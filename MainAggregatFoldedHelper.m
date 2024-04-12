function [real_porosity, K1, K2] = MainAggregatFoldedHelper( POR )

createCfile = false;

TestNum = 10;
% vecElem = (10 : 10 + TestNum);
% vecPDE  = zeros(TestNum+1,1);
% vecCAM  = zeros(TestNum+1,1);
% vecALL  = zeros(TestNum+1,1);

vecRun  = (1 : TestNum); 
vecEig  = zeros(2*TestNum,1);

format long
if createCfile
    oldMatrix = [1, 0; 0 0];
    fid = fopen('HomogenizedTensor.C', 'wt' );
    fprintf( fid , '/*\n * file: HomogenizedTensor.C\n');
    fprintf( fid , ' * author: Andreas Rupp\n');
    fprintf( fid , ' * Note: This file was automatically created using a MATLAB tool written by Andreas Rupp.\n');
    fprintf( fid , ' * For any information mail andreas.rupp@fau.de\n */\n\n\n');
    fprintf( fid , '#include "m++.h"\n\n\n');
    fprintf( fid , 'SmallMatrix homTensor(double t, Point p)  {\n');
end

for numberX = 1 : TestNum

% global bulkVector bioMvector

PDEtime = 0;
CAMtime = 0;

%% Definition of parameters                                                 %% Can be changed
NX          = 16;                       % Number of horicontal Elements
NZd         = NX;                       % Number of vertical Elements
porosity    = POR;                      % Porosity
numSmallPor = 0.5;                      % Approximate percentage of small pores
numBigPor   = 0.5;                      % Approximate percentage of big pores

tau         = 1;                        % End-Time
epsilon     = 1;                        % Everything < epsilon is defined as zero
numInnerIt  = 10;                       % Number of Diffusion/Absorption Iterations
numOuterIt  = 50;                       % Number of Transport/Growth Iterations

p           = 1;                        % Order of Polynomials for LDG Disc.
ord         = 4;                        % Order of Gaussian Integration Rule
eta         = 1;                        % Penalty Parameter for LDG

startConBio = 10;                       % Initial Concentration of Biomass
maxConcBio  = 100;                      % Max. Concentration of Biomass -> Growth
minConcBio  = 10;                       % Min. Concentration of Biomass -> Shrinkage

f_uptakeAgent = 100;                     % Growth Rate of Agent
f_uptakeAir = -epsilon/2;               % Amount of consumed Air when Agent grows
f_decay     = -0.08;                      % Degeneration Rate of Agent
f_bioMgrow  = 10;                       % Growth Rate of Biomass
f_bioMdecr  = -0.1;                       % Degeneration Rate of Biomass
f_bioUptAir = -epsilon/3;               % Amount of Air consumed, when Agent grows
f_uptakeBacteria = 20;
f_decayBac = -epsilon/10;
f_oxyBac   = -epsilon/2;

if createCfile
    fprintf( fid , '\tdouble numSteps = %d.;\n\n' , numOuterIt );
    fprintf( fid , '\tSmallMatrix K(3,3);\n' );
    fprintf( fid , '\t\t\t\t\t\t\t\t\tK[0][2] = 0.;\n' );
    fprintf( fid , '\t\t\t\t\t\t\t\t\tK[1][2] = 0.;\n' );
    fprintf( fid , '\tK[2][0] = 0.;\tK[2][1] = 0.;\tK[2][2] = 1.;\n\n' );
    fprintf( fid , '\tif ( t < 0. )  {\n' );
    fprintf( fid , '\t\tK[0][0] = 1.;\tK[0][1] = 0.;\n' );
    fprintf( fid , '\t\tK[1][0] = 0.;\tK[1][1] = 1.;\n' );
    fprintf( fid , '\t}' );
end

%% Fixed Parameters (of general framework)
NZu         = 0;
width       = NX;
intBound    = NZd * ones(NX+1,1);
upBound     = NZd * ones(NX+1,1);

zero        = @(x,y) x-x + 0;
one         = @(x,y) x-x + 1;
% uZero       = @(x,y) 100 * (NX/8 < y) .* (y < NX/2) .* (NX/3 < x) .* (x < 2*NX/3) + 5;
% uZero       = @(x,y) 10 * one(x,y);
% bacZero     = @(x,y) 10 * (x > 4) .* (x < 6) .* (y > 4) .* (y < 6);
% uZero       = @(x,y) 4 * x .* (NX - x);
uZero       = @(x,y) 100 / NX^2 * x .* (NX-x);
% bacZero     = zero;
% bacZero     = @(x,y) 20 * (x > 6) .* (x < 10) .* (y > 6) .* (y < 10);
bacZero     = @(x,y) 100 / NX^2 * x .* (NX-x);
% bacZero     = @(x,y) (1 - ((x < 35) .* (x > 29) .* (y < 40) .* (y > 24))) .* bacZero(x,y);
q1Zero      = zero;
q2Zero      = zero;
q1ZeroBac   = zero;
q2ZeroBac   = zero;
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
concAgent       = 10*ones( g.numE , 1 );    % Initial Concentration of Agent
NX = g.NX;

%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed
%bulkVector      = randi([0,1], g.numTsub, 1); % 1 is bulk, 0 is air     % Vector containing Bulk Distribution
% bulkVector = [ones(1,1); zeros(7,1); ones(10,1); zeros(7,1)];
% 
% createMeshOutput(bulkVector, 0);
% 
% bulkVector = [ zeros(NX^2/4, 1 ); ...
%     repmat( [zeros(7*NX/16,1); ones(NX/8,1); zeros(7*NX/16,1)], 3*NX/16, 1); ...
%     repmat( [zeros(3*NX/8,1); ones(5*NX/64,1); ones(3*NX/32,1); ones(5*NX/64,1); zeros(3*NX/8,1)], NX/8, 1); ...
%     repmat( [zeros(7*NX/16,1); ones(NX/8,1); zeros(7*NX/16,1)], 3*NX/16, 1); ...
%     zeros(NX^2/4, 1 ) ];
% bioMvector = zeros(size(bulkVector));

bulkVector = [ zeros(4*16,1); ...
    zeros(4,1); ones(8,1); zeros(4,1); ...
    zeros(4,1); ones(8,1); zeros(4,1); ...
    repmat( [zeros(4,1); ones(2,1); zeros(4,1); ones(2,1); zeros(4,1)], 6, 1 ); ...
    zeros(4*16,1) ];

% helperPor = porosity;
% while true
%     bulkVector = createBulkVector(g, floor(NX*NX*helperPor*numBigPor), floor(NX*NX*helperPor*numSmallPor));
%     if sum(~bulkVector) < NX^2 * (porosity - 0.0025)
%         helperPor = helperPor + 0.0025;
%     elseif sum(~bulkVector) > NX^2 * (porosity + 0.0025)
%         helperPor = helperPor - 0.0025;
%     else
%         break
%     end
% end
bioMvector = zeros(size(bulkVector));

%bulkVector = [zeros(44,1); 1; zeros(10,1); 1; zeros(10,1); 1; zeros(10,1); 1; zeros(22,1)];

real_porosity = sum(~bulkVector) / (NX^2)

% Note all data must be adapted manually in the following function (except
% for the bulkVector)

%bulkVector = [zeros(30,1); ones(40,1); zeros(30,1)];

% bulkVector = zeros(NX^2,1);
% bioMvector = zeros(size(bulkVector));

% bioMvector = [ zeros(3*NX,1); zeros(3,1); ones(4,1); zeros(2*3,1); 1; zeros(2,1); 1; zeros(2*3,1); 1; zeros(2,1); 1; zeros(2*3,1); ...
%     ones(4,1); zeros(3,1); zeros(3*NX,1) ];
% size(bulkVector)
% size(bioMvector)

% bulkVector = zeros(100,1);
% bioMvector = zeros(10^2,1);

% bioMvector      = randi([0,1], g.numTsub, 1); % 1 is bioM, 0 is air     % Vector containing Biomass Distribution
% bioMvector      = zeros(size(bulkVector));
% bioMvector      = bioMvector .* ~bulkVector;                            % Ensuring that Bulk and Biomass are disjoint

% bulkVector = zeros(g.numTsub, 1);
% for i = 1 : NZd
%     if mod(i,2) == 1
%         for j = 1 : NX
%             if mod(j,2) == 1
%                 bulkVector( (i-1)*NX + j ) = 1;
%             end  % if
%         end  % for j
%     end  % if
% end  % for
% bioMvector = zeros(size(bulkVector));

% global startBulkVector startBioMvector
% 
% if length(startBulkVector) == NX^2
%     bulkVector = startBulkVector;
%     bioMvector = startBioMvector;
% end  % if
% startBulkVector = bulkVector;
% startBioMvector = bioMvector;
airVector = ~bulkVector;
% porosity = sum(airVector) / g.numT;
% porosityBio = (sum(airVector) - sum(bioMvector)) / g.numT;

airVector = ~bulkVector;
bioMvec0airVec = bioMvector(airVector);
K1 = HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector);
if sum(sum(isnan(K1))) > 0
    K2 = 0;
    return;
end
% HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)))
% eigs(K)
numSolid = sum(bulkVector);


if createCfile && (sum(sum(isnan(K_hom))) == 0)
    fprintf( fid , '  else if ( t < %f )  {\n' , 1 / numOuterIt );
    printMatrix = K_hom;
    oldMatrix = K_hom;
    fprintf( fid , '\t\tK[0][0] = %f;\tK[0][1] = %f;\n' , printMatrix(1,1), printMatrix(1,2) );
    fprintf( fid , '\t\tK[1][0] = %f;\tK[1][1] = %f;\n\t}' , printMatrix(2,1), printMatrix(2,2) );
end    

%% Defining Initial Concentrations of Biomass
bioConcVector   = startConBio * bioMvector;
biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
% visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', 0);
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

for k = 1 : numOuterIt
    
    assert(numSolid == sum(bulkVector), 'solid Bloecke stimmen nich')
    
    %% Creating Vectors characterizing the Type of Edges
    idE             = zeros(g.numE, 1);
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
        bacLagr     = projectDG2LagrangeSub( UDGbac );
%         visualizeDataSub(g, uLagr, 'u', 'solu', 0);
%         visualizeDataSub(g, bacLagr, 'bac', 'solBac', 0);
    end  % if

    %% End of first PDE-time
    PDEtime = PDEtime + toc;

    %% Inner Loop for Diffusion and Growth/Shrinkage, ... (without Movement)
    
    for i = 1 : numInnerIt
        
        %% Starting PDE-time
        tic;
        
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
        
        %% End of PDE-time
        PDEtime = PDEtime + toc;
        
        %% Start of CAM-time
        tic;
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
        bioConcVector(bioMnoAir) = bioConcVector(bioMnoAir) + tau * f_bioMdecr * bioConcVector(bioMnoAir);
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
        %% End of CAM-time
        CAMtime = CAMtime + toc;
        
%        concAir = computeConcAir(g, uDG, ord);
%        concAir = concAir < epsilon;
%        uDG(concAir,:) = 0;
    end  % for i
% norm(concAgent,inf)
    %% Start of CAM-time
    tic;
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
    helper1 = [ candidates - NX , candidates - 1 , candidates , candidates + 1 , candidates + NX ];
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
    
    
    candidates = find(OverHead);
    helper1 = [  candidates - 2 * NX , candidates - NX - 1 , candidates - NX , candidates - NX + 1 , ...
            candidates - 2 , candidates - 1 , candidates , candidates + 1 , candidates + 2 , ...
            candidates + NX - 1 , candidates + NX , candidates + NX + 1 , candidates + 2 * NX ];
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
    concAgent = concAgent .* (concAgent > epsilon);
    NoBioMassVec = ~bioMvector;
    
    bulk0T = zeros( g.numT , 4 );
    bulk0T( g.NX + 1 : end , 1 ) = bulkVector( 1 : end-g.NX );
    bulk0T( 1 : end - g.NX , 2 ) = bulkVector( g.NX+1 : end );
    bulk0T( 1 : end - 1 , 3 ) = bulkVector( 2 : end );
    bulk0T( 2 : end , 4 ) = bulkVector( 1 : end - 1 );
    
    candidates = bulkVector .* (concAgent(g.E0T(:,1)) == 0 | ~bulk0T(:,1)) .* (concAgent(g.E0T(:,2)) == 0 | ~bulk0T(:,2)) ...
                    .* (concAgent(g.E0T(:,3)) == 0 | ~bulk0T(:,3)) .* (concAgent(g.E0T(:,4)) == 0 | ~bulk0T(:,4)) .* (bioMvector == 0);
    candidates(2:end)       = candidates(2:end) .* NoBioMassVec(1:end-1);
    candidates(NX+1:end)    = candidates(NX+1:end) .* NoBioMassVec(1:end-NX);
    candidates(1:end-NX)    = candidates(1:end-NX) .* NoBioMassVec(NX+1:end);
    candidates(1:end-1)     = candidates(1:end-1) .* NoBioMassVec(2:end);
                
    candidates = find(candidates);
    
    helper1 = [  candidates - 2 * NX , candidates - NX - 1 , candidates - NX , candidates - NX + 1 , ...
                candidates - 2 , candidates - 1 , candidates , candidates + 1 , candidates + 2 , ...
                candidates + NX - 1 , candidates + NX , candidates + NX + 1 , candidates + 2 * NX ];
    helper2 = (helper1 > 0) .* (helper1 < g.numT + 1);
    helper2 = helper2 > 0;
    aims = zeros(length(candidates),13);
    aims(helper2) = helper1(helper2);
    aims = reshape(aims, length(candidates), 13)';
    values = zeros(size(aims));
    
    for i = 1 : size(aims,1)
        for j = 1 : size(aims,2)
            if aims(i,j) > 0 && i ~= 7
                aims(i,j) = aims(i,j) .* airVector(aims(i,j)) .* ~bioMvector(aims(i,j));
            end
            if aims(i,j) > 0
                values(i,j) = concAgent(g.E0T(aims(i,j),1)) + concAgent(g.E0T(aims(i,j),2)) ...
                    + concAgent(g.E0T(aims(i,j),3)) + concAgent(g.E0T(aims(i,j),4)) ...
                    + idEbulk(g.E0T(aims(i,j),1)) + idEbulk(g.E0T(aims(i,j),2)) ...
                    + idEbulk(g.E0T(aims(i,j),3)) + idEbulk(g.E0T(aims(i,j),4));
            end
            if i == 7
                values(i,j) = concAgent(g.E0T(aims(i,j),1)) * idEbulk2(g.E0T(aims(i,j),1)) + concAgent(g.E0T(aims(i,j),2)) * idEbulk2(g.E0T(aims(i,j),2)) ...
                    + concAgent(g.E0T(aims(i,j),3)) * idEbulk2(g.E0T(aims(i,j),3)) + concAgent(g.E0T(aims(i,j),4)) * idEbulk2(g.E0T(aims(i,j),4)) ...
                    + idEbulk2(g.E0T(aims(i,j),1)) + idEbulk2(g.E0T(aims(i,j),2)) ...
                    + idEbulk2(g.E0T(aims(i,j),3)) + idEbulk2(g.E0T(aims(i,j),4));
            end
            if (i == 3) || (i == 6) || (i == 8) || (i == 11)
                values(i,j) = values(i,j) - 1;
                values(i,j) = max(values(i,j), 0);
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

%    aims = max(aims)';
    [n, bin] = histc(aims, unique(aims));
    multiple = find(n > 1);
    index    = find(ismember(bin, multiple));
    
    for i = 1 : length(index)-1
        if aims(index(i)) ~= 0
%             for j = i+1 : length(index)
%                 if aims(index(i)) == aims(index(j))
%                     candidates(index(j)) = 0;
%                     aims(index(j)) = 0;
%                 end  % if equal
%             end  % for j
            aims(index(i+1:end)) = aims(index(i+1:end)) .* (aims(index(i+1:end)) ~= aims(index(i)));
        end  % if zero
    end  % for i
    
%    candidates(index(2:end)) = 0;
%    aims(index(2:end)) = 0;
    candidates = candidates(aims ~= 0);
    aims = aims(aims ~= 0);
    
    assert(length(candidates) == length(aims), 'Laenge passt nicht');
    
    neighbours = [aims-NX, aims - 1, aims + 1, aims + NX];
    conflicts = zeros(length(aims), 4);
    for n = 1 : 4
        neighbours(:,n) = neighbours(:,n) .* (neighbours(:,n) < g.numT + 1);
        for i = 1 : length(aims)
            if neighbours(i,n) > 0 && bulkVector(neighbours(i,n)) == 1 && candidates(i) ~= neighbours(i,n) && any(candidates == neighbours(i,n))
                conflicts(i,n) = 1;                
            end
        end 
    end

   while sum(sum(conflicts)) ~= 0
       problems = neighbours(conflicts == 1);
       index = randi(length(problems));
       mf = problems(index);
       cand = find(candidates == mf);
       aims(cand) = 0;
       candidates(cand) = 0;
       conflicts(cand,:) = [0,0,0,0];
       cand = find(neighbours == mf);
       [a,b] = size(conflicts);
       conflicts = reshape(conflicts, a*b, 1);
       conflicts(cand) = 0;
       conflicts = reshape(conflicts, a, b);
   end
    
    candidates = candidates(aims ~= 0);
    aims = aims(aims ~= 0);
    
    bulkVector(candidates) = 0;
    bulkVector(aims) = 1;
    uDG(candidates,:) = uDG(aims, :);
    uDG(aims,:) = zeros(length(aims), N);
    q1DG(candidates,:) = q1DG(aims, :);
    q1DG(aims,:) = zeros(length(aims), N);
    q2DG(candidates,:) = q2DG(aims, :);
    q2DG(aims,:) = zeros(length(aims), N);
    uDGbac(candidates,:) = uDGbac(aims, :);
    uDGbac(aims,:) = zeros(length(aims), N);
    q1DGbac(candidates,:) = q1DGbac(aims, :);
    q1DGbac(aims,:) = zeros(length(aims), N);
    q2DGbac(candidates,:) = q2DGbac(aims, :);
    q2DGbac(aims,:) = zeros(length(aims), N);
    
%     if k == 1 | k == numOuterIt
%         numAir
%     end

    candidates = bulkVector .* ( (concAgent(g.E0T(:,1)) == 0 | ~bulk0T(:,1)) + (concAgent(g.E0T(:,2)) == 0 | ~bulk0T(:,2)) ...
                    + (concAgent(g.E0T(:,3)) == 0 | ~bulk0T(:,3)) + (concAgent(g.E0T(:,4)) == 0 | ~bulk0T(:,4)) );
    candidates = (candidates == 3);
    
    markE0TE0T = cell(4,1);
    for i = 1 : 4
        markE0TE0T{i} = sparse( g.numT , g.numT );
    end  % for
    markE0TE0T{1}( sub2ind([g.numT,g.numT], NX+1:g.numT , 1:g.numT-NX) ) = 1;
    markE0TE0T{1} = bsxfun(@times, (concAgent(g.E0T(:,1)) > 0), markE0TE0T{1});
    markE0TE0T{2}( sub2ind([g.numT,g.numT], 1:g.numT-NX , NX+1:g.numT) ) = 1;
    markE0TE0T{2} = bsxfun(@times, (concAgent(g.E0T(:,2)) > 0), markE0TE0T{2});
    markE0TE0T{3}( sub2ind([g.numT,g.numT], 1:g.numT-1 , 2:g.numT) ) = 1;
    markE0TE0T{3}( sub2ind([g.numT,g.numT], NX:NX:g.numT-NX , NX+1:NX:g.numT-NX+1) ) = 0;
    markE0TE0T{3} = bsxfun(@times, (concAgent(g.E0T(:,3)) > 0), markE0TE0T{3});
    markE0TE0T{4}( sub2ind([g.numT,g.numT], 2:g.numT , 1:g.numT-1) ) = 1;
    markE0TE0T{4}( sub2ind([g.numT,g.numT], NX+1:NX:g.numT-NX+1 , NX:NX:g.numT-NX) ) = 0;
    markE0TE0T{4} = bsxfun(@times, (concAgent(g.E0T(:,4)) > 0), markE0TE0T{4});
    E0TE0T = markE0TE0T{1} + markE0TE0T{2} + markE0TE0T{3} + markE0TE0T{4};
    E0TE0T = E0TE0T(candidates,candidates);
    [cands1, cands2] = find(E0TE0T);
%     cands = find(candidates);
 
    if ~isempty(cands1)
       for i = 1 : length(cands1)
           for j = 1 : length(cands1)
               if cands1(i) == cands2(j)
                cands2(j) = 0;
               end
           end
       end
    end
    
    cands1 = cands1(cands2 > 0);
    cands2 = cands2(cands2 > 0);
 
    helper1 = [  cands1 - 2 * NX , cands1 - NX - 1 , cands1 - NX , cands1 - NX + 1 , ...
                cands1 - 2 , cands1 - 1 , cands1 , cands1 + 1 , cands1 + 2 , ...
                cands1 + NX - 1 , cands1 + NX , cands1 + NX + 1 , cands1 + 2 * NX ];
    helper1 = (helper1 > 0) .* (helper1 < g.numT + 1);
    helper1 = helper1 > 0;
    
    aims    = zeros(size(helper1,1),2);
    
    if ~isempty(helper1)
       for i = 1 : size(helper1,1)
           for j = 1 : size(helper1,2)
               helper1(i,j) = helper1(i,j) * airVector(helper1(i,j)) * ~bioMvector(helper1(i,j));
               if helper(i,j) ~= 0
                   if helper(i,j)-NX > 0 && airVector(helper1(i,j)-NX) && ~bioMvector(helper1(i,j)-NX)
                       aims(i,1) = helper1(i,j);
                       aims(i,2) = helper1(i,j)-NX;
                   elseif helper(i,j)-1 > 0 && airVector(helper1(i,j)-1) && ~bioMvector(helper1(i,j)-1)
                       aims(i,1) = helper1(i,j);
                       aims(i,2) = helper1(i,j)-1;
                   elseif helper(i,j)+1 < g.numT && airVector(helper1(i,j)+1) && ~bioMvector(helper1(i,j)+1)
                       aims(i,1) = helper1(i,j);
                       aims(i,2) = helper1(i,j)+1;
                   elseif helper(i,j)+NX < g.numT && airVector(helper1(i,j)+NX) && ~bioMvector(helper1(i,j)+NX)
                       aims(i,1) = helper1(i,j);
                       aims(i,2) = helper1(i,j)+NX;    
                   end
               end
           end
       end
    end
    
    aims(:,1) = aims(:,1) .* ~bioMvector(aims(:,1)) .* ~bioMvector(aims(:,2));

    cands1 = cands1(aims(:,1) > 0);
    cands2 = cands2(aims(:,1) > 0);
    aims   = aims(aims(:,1) > 0,:);

    bulkVector(cands1) = 0;
    bulkVector(cands2) = 0;
    bulkVector(aims(:,1)) = 1;
    bulkVector(aims(:,2)) = 1;    
    uDG(cands1,:) = uDG(aims(:,1), :);
    uDG(cands2,:) = uDG(aims(:,2), :);
    uDG(aims(:,1),:) = zeros(length(aims(:,1)), N);
    uDG(aims(:,2),:) = zeros(length(aims(:,2)), N);
    q1DG(cands1,:) = q1DG(aims(:,1), :);
    q1DG(cands2,:) = q1DG(aims(:,2), :);
    q1DG(aims(:,1),:) = zeros(length(aims(:,1)), N);
    q1DG(aims(:,1),:) = zeros(length(aims(:,2)), N);
    q2DG(cands1,:) = q2DG(aims(:,1), :);
    q2DG(cands2,:) = q2DG(aims(:,2), :);
    q2DG(aims(:,1),:) = zeros(length(aims(:,1)), N);
    q2DG(aims(:,1),:) = zeros(length(aims(:,2)), N);
    uDGbac(cands1,:) = uDGbac(aims(:,1), :);
    uDGbac(cands2,:) = uDGbac(aims(:,2), :);
    uDGbac(aims(:,1),:) = zeros(length(aims(:,1)), N);
    uDGbac(aims(:,2),:) = zeros(length(aims(:,2)), N);
    q1DGbac(cands1,:) = q1DGbac(aims(:,1), :);
    q1DGbac(cands2,:) = q1DGbac(aims(:,2), :);
    q1DGbac(aims(:,1),:) = zeros(length(aims(:,1)), N);
    q1DGbac(aims(:,1),:) = zeros(length(aims(:,2)), N);
    q2DGbac(cands1,:) = q2DGbac(aims(:,1), :);
    q2DGbac(cands2,:) = q2DGbac(aims(:,2), :);
    q2DGbac(aims(:,1),:) = zeros(length(aims(:,1)), N);
    q2DGbac(aims(:,1),:) = zeros(length(aims(:,2)), N);
    
    %% Pairs of Bulk
    
    numBulkOld = sum(bulkVector);
    
    bulk0T( g.NX + 1 : end , 1 ) = bulkVector( 1 : end-g.NX );
    bulk0T( 1 : end - g.NX , 2 ) = bulkVector( g.NX+1 : end );
    bulk0T( 1 : end - 1 , 3 ) = bulkVector( 2 : end );
    bulk0T( 2 : end , 4 ) = bulkVector( 1 : end - 1 );
    
    idE = zeros(g.numE,1);
    
    idE(g.E0T(bulkVector == 1,1)) = idE(g.E0T(bulkVector == 1,1)) + 1;
    idE(g.E0T(bulkVector == 1,2)) = idE(g.E0T(bulkVector == 1,2)) + 1;
    idE(g.E0T(bulkVector == 1,3)) = idE(g.E0T(bulkVector == 1,3)) + 1;
    idE(g.E0T(bulkVector == 1,4)) = idE(g.E0T(bulkVector == 1,4)) + 1;
    idE = idE + markEbdr;

    idEairy         = idE == 0;
    idEbulk         = ~idEairy;
    idEbulk2        = idE == 2;
    
    candidates = bulkVector .* (concAgent(g.E0T(:,1)) ~= 0 & bulk0T(:,1) & (concAgent(g.E0T(:,2)) == 0 | ~bulk0T(:,2)) & (concAgent(g.E0T(:,3)) == 0 | ~bulk0T(:,3)) & (concAgent(g.E0T(:,4)) == 0 | ~bulk0T(:,4)) ) ...
        + bulkVector .* (concAgent(g.E0T(:,4)) ~= 0 & bulk0T(:,4) & (concAgent(g.E0T(:,1)) == 0 | ~bulk0T(:,1)) & (concAgent(g.E0T(:,2)) == 0 | ~bulk0T(:,2)) & (concAgent(g.E0T(:,3)) == 0 | ~bulk0T(:,3)) );
    candidates = find(candidates);
    
    concAgent(idEairy) = 0;
    
    if ~isempty(candidates)
        for i = 1 : length(candidates)
            if bulk0T(candidates(i),1)
                if candidates(i)-NX <1
                    candidates(i) = 0;
                else
                    candidates(i) = candidates(i) * ( concAgent(g.E0T(candidates(i) - NX,1)) == 0 | ~bulk0T(candidates(i) - NX,1) ) ...
                        * ( concAgent(g.E0T(candidates(i) - NX,3)) == 0 | ~bulk0T(candidates(i) - NX,3) ) ...
                        * ( concAgent(g.E0T(candidates(i) - NX,4)) == 0 | ~bulk0T(candidates(i) - NX,4) );
                end
            else
                if candidates(i) < 2
                    candidates(i) = 0;
                else
                    candidates(i) = candidates(i) * ( concAgent(g.E0T(candidates(i) - 1,1)) == 0 | ~bulk0T(candidates(i) - 1,1) ) ...
                        * ( concAgent(g.E0T(candidates(i) - 1,2)) == 0 | ~bulk0T(candidates(i) - 1,2) ) ...
                        * ( concAgent(g.E0T(candidates(i) - 1,4)) == 0 | ~bulk0T(candidates(i) - 1,4) );
                end
            end  % if bulk0T
        end  % for i
    end  % if ~isempty(candidates)
    
    candidates = candidates(candidates > 0);
    
    if ~isempty(candidates)
        helping = [ candidates - 2*NX , candidates - NX, candidates - 2, candidates - 1, candidates, ...
                    candidates + 1, candidates + 2, candidates + NX, candidates + 2 * NX ];
        helping = helping .* (helping > 0);
        helping = helping .* (helping < g.numT + 1);
        aims = zeros(size(helping));
        for i = 1 : length(candidates)
            if bulk0T(candidates(i),1)
                for j = [1, 3, 4, 6, 7, 9]
                    if helping(i,j) - g.NX < 1
                        helping(i,j) = 0;
                    elseif helping(i,j) ~= 0
                        helping(i,j) = helping(i,j) * ~bulkVector(helping(i,j)) * ~bulkVector(helping(i,j) - g.NX);
                    end  % if
                end  % for j
                if helping(i,2) - g.NX < 1
                    helping(i,2) = 0;
                elseif helping(i,2) ~= 0
                    helping(i,2) = helping(i,2) * ~bulkVector(helping(i,2) - g.NX);
                end
                if helping(i,8) + 1 > g.numT
                    helping(i,8) = 0;
                elseif helping(i,8) ~= 0
                    helping(i,8) = helping(i,8) * ~bulkVector(helping(i,8));
                end
                for j = [1, 2, 3, 4, 6, 7, 8, 9]
                    if helping(i,j) ~= 0
                        aims(i,j) = concAgent(g.E0T(helping(i,j),2))  * idEbulk(g.E0T(helping(i,j),2))...
                            + concAgent(g.E0T(helping(i,j),3)) * idEbulk(g.E0T(helping(i,j),3)) + concAgent(g.E0T(helping(i,j),4)) * idEbulk(g.E0T(helping(i,j),4)) ...
                            + idEbulk(g.E0T(helping(i,j),2)) ...
                            + idEbulk(g.E0T(helping(i,j),3)) + idEbulk(g.E0T(helping(i,j),4)) ...
                            + concAgent(g.E0T(helping(i,j)-g.NX,1)) * idEbulk(g.E0T(helping(i,j)-g.NX,1)) ...
                            + concAgent(g.E0T(helping(i,j)-g.NX,3)) * idEbulk(g.E0T(helping(i,j)-g.NX,3)) + concAgent(g.E0T(helping(i,j)-g.NX,4)) * idEbulk(g.E0T(helping(i,j)-g.NX,4)) ...
                            + idEbulk(g.E0T(helping(i,j)-g.NX,1)) ...
                            + idEbulk(g.E0T(helping(i,j)-g.NX,3)) + idEbulk(g.E0T(helping(i,j)-g.NX,4));
                    end
                end  % for j
                if helping(i,5) ~= 0
                    aims(i,5) = concAgent(g.E0T(helping(i,5),2)) * idEbulk2(g.E0T(helping(i,5),2)) ...
                        + concAgent(g.E0T(helping(i,5),3)) * idEbulk2(g.E0T(helping(i,5),3)) + concAgent(g.E0T(helping(i,5),4)) * idEbulk2(g.E0T(helping(i,5),4)) ...
                        + idEbulk2(g.E0T(helping(i,5),2)) ...
                        + idEbulk2(g.E0T(helping(i,5),3)) + idEbulk2(g.E0T(helping(i,5),4)) ...
                        + concAgent(g.E0T(helping(i,5)-g.NX,1)) * idEbulk2(g.E0T(helping(i,5)-g.NX,1)) ...
                        + concAgent(g.E0T(helping(i,5)-g.NX,3)) * idEbulk2(g.E0T(helping(i,5)-g.NX,3)) + concAgent(g.E0T(helping(i,5)-g.NX,4)) * idEbulk2(g.E0T(helping(i,5)-g.NX,4)) ...
                        + idEbulk2(g.E0T(helping(i,5)-g.NX,1)) ...
                        + idEbulk2(g.E0T(helping(i,5)-g.NX,3)) + idEbulk2(g.E0T(helping(i,5)-g.NX,4));
                end  % if
            else
                for j = [1, 2, 3, 7, 8, 9]
                    if helping(i,j) -1 < 1
                        helping(i,j) = 0;
                    elseif helping(i,j) ~= 0
                        helping(i,j) = helping(i,j) * ~bulkVector(helping(i,j)) * ~bulkVector(helping(i,j) - 1);
                    end  % if
                end  % for j
                if helping(i,4) - 1 < 1
                    helping(i,4) = 0;
                elseif helping(i,j) ~= 0
                    helping(i,4) = helping(i,4) * ~bulkVector(helping(i,4) - 1);
                end
                if helping(i,6) - 1 < 1
                    helping(i,6) = 0;
                elseif helping(i,6) ~= 0
                    helping(i,6) = helping(i,6) * ~bulkVector(helping(i,6));
                end
                for j = [1, 2, 3, 4, 6, 7, 8, 9]
                    if helping(i,j) ~= 0
                        aims(i,j) = concAgent(g.E0T(helping(i,j),1)) * idEbulk(g.E0T(helping(i,j),1)) + concAgent(g.E0T(helping(i,j),2)) * idEbulk(g.E0T(helping(i,j),2))...
                            + concAgent(g.E0T(helping(i,j),3)) * idEbulk(g.E0T(helping(i,j),3)) ...
                            + idEbulk(g.E0T(helping(i,j),1)) + idEbulk(g.E0T(helping(i,j),2)) ...
                            + idEbulk(g.E0T(helping(i,j),3)) ...
                            + concAgent(g.E0T(helping(i,j)-1,1)) * idEbulk(g.E0T(helping(i,j)-1,1)) + concAgent(g.E0T(helping(i,j)-1,2)) * idEbulk(g.E0T(helping(i,j)-1,2)) ...
                            + concAgent(g.E0T(helping(i,j)-1,4)) * idEbulk(g.E0T(helping(i,j)-1,4)) ...
                            + idEbulk(g.E0T(helping(i,j)-1,1)) + idEbulk(g.E0T(helping(i,j)-1,2)) ...
                            + idEbulk(g.E0T(helping(i,j)-1,4));
                    end
                end  % for j
                if helping(i,5) ~= 0
                    aims(i,5) = concAgent(g.E0T(helping(i,5),1)) * idEbulk2(g.E0T(helping(i,5),1)) + concAgent(g.E0T(helping(i,5),2)) * idEbulk2(g.E0T(helping(i,5),2)) ...
                        + concAgent(g.E0T(helping(i,5),3)) * idEbulk2(g.E0T(helping(i,5),3))...
                        + idEbulk2(g.E0T(helping(i,5),1)) + idEbulk2(g.E0T(helping(i,5),2)) ...
                        + idEbulk2(g.E0T(helping(i,5),3)) ...
                        + concAgent(g.E0T(helping(i,5)-1,1)) * idEbulk2(g.E0T(helping(i,5)-1,1)) + concAgent(g.E0T(helping(i,5)-1,2)) * idEbulk2(g.E0T(helping(i,5)-1,2))...
                        + concAgent(g.E0T(helping(i,5)-1,4)) * idEbulk2(g.E0T(helping(i,5)-1,4))...
                        + idEbulk2(g.E0T(helping(i,5)-1,1)) + idEbulk2(g.E0T(helping(i,5)-1,2)) ...
                        + idEbulk2(g.E0T(helping(i,5)-1,4));
                end  % if
            end  % if bulk0T
        end  % for i
    end  % if ~ismepty candidates
    
    
    
    [maximo, indMax] = max(aims');
    aims2 = zeros(size(indMax'));
    for i = 1 : length(indMax)
%         aims2(i) = aims(indMax(i),i);
      indices = find( maximo(i) == aims(i,:) );
      chooser = randi(length(indices));
      if length(indices) == 1 && indMax(i) ~= indices(chooser)
          error('IDIOT')
      end
      aims2(i) = helping(i, indices(chooser));
    end

    aims = aims2;
    
    candidates = candidates(aims > 0);
    aims = aims(aims > 0);
    
    uDGold = uDG;
    bulkVectorOld = bulkVector;
    
    for i = 1 : length(aims)
        bulkVector(candidates(i)) = 0;
        uDG(candidates(i),:) = uDG(aims(i),:);
        if bulk0T(candidates(i),1)
            bulkVector(candidates(i) - g.NX ) = 0;
            bulkVector(aims(i) - g.NX) = 1;
            uDG(candidates(i) - g.NX,:) = uDG(aims(i)-g.NX,:);
             if uDG(candidates(i),:) == zeros(1,N)
               uDG(candidates(i),:) = uDG(aims(i) - g.NX,:); 
            end
            uDG(aims(i) - g.NX, :) = zeros(1,N);
        else
            bulkVector(candidates(i) - 1 ) = 0;
            bulkVector(aims(i) - 1) = 1;
            uDG(candidates(i) - 1,:) = uDG(aims(i) - 1,:);
            uDG(aims(i) - 1, :) = zeros(1,N);
            if uDG(candidates(i),:) == zeros(1,N)
               uDG(candidates(i),:) = uDG(aims(i),:); 
            end
        end  % if
        uDG(aims(i),:) = zeros(1, N);
        bulkVector(aims(i)) = 1;
    end  % for i

    if (numBulkOld ~= sum(bulkVector) || sum(bioMvector .* bulkVector ~= 0))
        uDG = uDGold;
        bulkVector = bulkVectorOld;
        i = 0;
        while (numBulkOld ~= sum(bulkVector) || sum(bioMvector .* bulkVector ~= 0)) 
            i = i + 1;
            bulkVector(candidates(i)) = 0;
            uDG(candidates(i),:) = uDG(aims(i),:);
            if bulk0T(candidates(i),1)
                bulkVector(candidates(i) - g.NX ) = 0;
                bulkVector(aims(i) - g.NX) = 1;
                uDG(candidates(i) - g.NX,:) = uDG(aims(i)-g.NX,:);
                if uDG(candidates(i),:) == zeros(1,N)
                    uDG(candidates(i),:) = uDG(aims(i) - g.NX,:); 
                end
                uDG(aims(i) - g.NX, :) = zeros(1,N);
            else
                bulkVector(candidates(i) - 1 ) = 0;
                bulkVector(aims(i) - 1) = 1;
                uDG(candidates(i) - 1,:) = uDG(aims(i) - 1,:);
                uDG(aims(i) - 1, :) = zeros(1,N);
                if uDG(candidates(i),:) == zeros(1,N)
                    uDG(candidates(i),:) = uDG(aims(i),:); 
                end
            end  % if
            uDG(aims(i),:) = zeros(1, N);
            bulkVector(aims(i)) = 1;
        end
        corrects = i - 1;
        uDG = uDGold;
        bulkVector = bulkVectorOld;
        for i = 1 : corrects
            bulkVector(candidates(i)) = 0;
            uDG(candidates(i),:) = uDG(aims(i),:);
            if bulk0T(candidates(i),1)
                bulkVector(candidates(i) - g.NX ) = 0;
                bulkVector(aims(i) - g.NX) = 1;
                uDG(candidates(i) - g.NX,:) = uDG(aims(i)-g.NX,:);
                if uDG(candidates(i),:) == zeros(1,N)
                    uDG(candidates(i),:) = uDG(aims(i) - g.NX,:); 
                end
                uDG(aims(i) - g.NX, :) = zeros(1,N);
            else
                bulkVector(candidates(i) - 1 ) = 0;
                bulkVector(aims(i) - 1) = 1;
                uDG(candidates(i) - 1,:) = uDG(aims(i) - 1,:);
                uDG(aims(i) - 1, :) = zeros(1,N);
                if uDG(candidates(i),:) == zeros(1,N)
                    uDG(candidates(i),:) = uDG(aims(i),:); 
                end
            end  % if
            uDG(aims(i),:) = zeros(1, N);
            bulkVector(aims(i)) = 1;
        end  % for i
    end  % if numBulkOld
    
    %% End of CAM-time
    CAMtime = CAMtime + toc;
    
%     porosity = sum(airVector) / g.numT;
%     porosityBio = (sum(airVector) - sum(bioMvector)) / g.numT;

    
    airVector = ~bulkVector;
    bioMvec0airVec = bioMvector(airVector);
%     K_hom = HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector)
%     HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)))
  
%     if i == 1 || i == numInnerIt
        uLagr       = projectDG2LagrangeSub( uDG );
        bacLagr     = projectDG2LagrangeSub( uDGbac );
%         visualizeDataSub(g, uLagr, 'u', 'solu', k);
%         visualizeDataSub(g, bacLagr, 'bac', 'solBac', k);
        biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
%         visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', k);
%     end

    if sum(bioMvector .* bulkVector ~= 0)
        error('NEEEEEEEEEEEIN!')
    end
    
    if createCfile
        if k == numOuterIt
            fprintf( fid , '  else  {\n' );
            if sum(sum(isnan(K_hom))) > 0
                printMatrix = oldMatrix;
            else
                printMatrix = K_hom;
            end
            fprintf( fid , '\t\tK[0][0] = %f;\tK[0][1] = %f;\n' , printMatrix(1,1), printMatrix(1,2) );
            fprintf( fid , '\t\tK[1][0] = %f;\tK[1][1] = %f;\n\t}' , printMatrix(2,1), printMatrix(2,2) );
        else
            if sum(sum(isnan(K_hom))) == 0
                fprintf( fid , '  else if ( t < %f )  {\n' , (k+1) / numOuterIt );
                printMatrix = K_hom;
                oldMatrix = K_hom;
                fprintf( fid , '\t\tK[0][0] = %f;\tK[0][1] = %f;\n' , printMatrix(1,1), printMatrix(1,2) );
                fprintf( fid , '\t\tK[1][0] = %f;\tK[1][1] = %f;\n\t}' , printMatrix(2,1), printMatrix(2,2) );
            end
        end
    end 
    
end  % for k

%     vecPDE(numberX-9) = PDEtime;
%     vecCAM(numberX-9) = CAMtime;
%     vecALL(numberX-9) = PDEtime + CAMtime;

% HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector);
% K = HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)));
% if sum(sum(isnan(K))) > 0
%     continue
% end
% a = eigs(K);
% 
% vecEig(2*numberX-1) = a(1);
% vecEig(2*numberX) = a(2);

end  % for numberX

if createCfile
    fprintf( fid , '\n\n\treturn K;\n}' );
    fclose( fid );
end

% figure('Name','PDE Solving Time')
% plot(vecElem, vecPDE)
% figure('Name','CAM Solving Time')
% plot(vecElem, vecCAM)
% figure('Name','Complete Time')
% plot(vecElem, vecALL)
% 
% ordPDE = polyfit(log(vecElem'), log(vecPDE), 1);
% ordCAM = polyfit(log(vecElem'), log(vecCAM), 1);
% ordALL = polyfit(log(vecElem'), log(vecALL), 1);
% ordPDE = ordPDE(1)
% ordCAM = ordCAM(1)
% ordALL = ordALL(1)

% airVector = ~bulkVector;
% bioMvec0airVec = bioMvector(airVector);
% 
% porosity = sum(airVector) / g.numT;
% porosityBio = (sum(airVector) - sum(bioMvector)) / g.numT;

% HomogenizationFunc(bulkVector, bioMvec0airVec, bioMvector)
K2 = HomogenizationFunc(bulkVector, zeros(size(bioMvec0airVec)), zeros(size(bulkVector)));
% eigs(K)

% plot(vecRun, vecEig(1:2:end), 'o', vecRun, vecEig(2:2:end), 'x')

% toc

end