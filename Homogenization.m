function Homogenization(  )

tic;

%% Definition of parameters                                                 %% Can be changed
NX          = 64;                       % Number of horicontal Elements
NZd         = 64;                       % Number of vertical Elements

tau         = 1;                        % End-Time
epsilon     = 1;                        % Everything < epsilon is defined as zero

p           = 1;                        % Order of Polynomials for LDG Disc.
ord         = 4;                        % Order of Gaussian Integration Rule
eta         = 1;                        % Penalty Parameter for LDG

%% Fixed Parameters (of general framework)
NZu         = 0;
width       = 1;
intBound    = 1 * ones(NX+1,1);
upBound     = 1 * ones(NX+1,1);

zero        = @(x,y) x-x + 0;
one         = @(x,y) x-x + 1;
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
% bulkVector      = randi([0,1], g.numTsub, 1); % 1 is bulk, 0 is air     % Vector containing Bulk Distribution

bulkVector = [ zeros(8*64,1); ...
    zeros(8,1); ones(5,1); zeros(51,1); ...
    zeros(8,1); ones(6,1); zeros(50,1); ...
    zeros(8,1); ones(7,1); zeros(49,1); ...
    zeros(8,1); ones(8,1); zeros(48,1); ...
    repmat([zeros(8,1);ones(9,1);zeros(48,1)], 39, 1); ...
    zeros(8,1); ones(9,1); zeros(8,1); ...
    zeros(48,1); ones(8,1); zeros(8,1); ...
    zeros(49,1); ones(7,1); zeros(8,1); ...
    zeros(50,1); ones(6,1); zeros(8,1); ...
    zeros(51,1); ones(5,1); zeros(8,1); ...
    zeros(8*64,1) ];
    
    length(bulkVector)
% bioMvector = zeros(100,1);

% bulkVector = zeros(64,1);
% bioMvector = zeros(64,1);

% bioMvector      = randi([0,1], g.numTsub, 1); % 1 is bioM, 0 is air     % Vector containing Biomass Distribution
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

%% Defining Initial Concentrations of Biomass
% bioConcVector   = startConBio * bioMvector;
% biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
% visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', 0);
% bioMvector      = bioMvector > 0;

%% Computation of Local Matrices for LDG Scheme
computeBasesOnQuad(p, ord);

[hatMc, hatMx]          = computeHatM(p, ord);
[hatGc, hatGx, hatGy]   = computeHatG(p, ord);
[hatHc, hatHx, hatHy]   = computeHatH(p, ord);
hatInt                  = computeHatInt(p, ord);
hatRdiag                = computeHatRdiag(p, ord);
hatRoffdiag             = computeHatRoffdiag(p, ord);
hatSdiag                = computeHatSdiag(p, ord);
hatSoffdiag             = computeHatSoffdiag(p, ord);

%% Creating Vectors characterizing the Type of Edges
idE             = zeros(g.numE, 1);
idE(g.E0T(bulkVector == 1,1)) = idE(g.E0T(bulkVector == 1,1)) + 1;
idE(g.E0T(bulkVector == 1,2)) = idE(g.E0T(bulkVector == 1,2)) + 1;
idE(g.E0T(bulkVector == 1,3)) = idE(g.E0T(bulkVector == 1,3)) + 1;
idE(g.E0T(bulkVector == 1,4)) = idE(g.E0T(bulkVector == 1,4)) + 1;
idE = idE + markEbdr;

idEneum         = idE == 1;
idEairy         = idE == 0;
% idEbulk         = ~(idEneum | idEairy);
markE0Tneum     = idEneum(g.E0T);
markE0Tairy     = idEairy(g.E0T);
airVector   = ~bulkVector;

markE0TneumSub = markE0Tneum(1:g.numTsub,:);
markE0TairySub = markE0Tairy(1:g.numTsub,:);
airVectorSub = airVector(1:g.numTsub);
markE0TneumSub = markE0TneumSub(airVectorSub,:);
markE0TairySub = markE0TairySub(airVectorSub,:);
numAir = norm(airVectorSub.*ones(size(airVectorSub)),1);

%% Assembling global Matrices (independent of u, q, concentrations)
globM       = assembleGlobMsub(g, hatMc, hatMx, airVectorSub, numAir);
globH       = assembleGlobHsub(g, hatHc, hatHx, hatHy, airVectorSub, numAir);
globQ       = assembleGlobQsub(g, markE0TairySub, hatSdiag, hatSoffdiag, airVectorSub, numAir);
globQN      = assembleGlobQNsub(g, markE0TneumSub, hatSdiag, airVectorSub, numAir);
globS       = assembleGlobSsub(g, markE0TairySub, hatSdiag, hatSoffdiag, eta, airVectorSub, numAir);

fDG         = projectAlg2DGsub(g, f, p, ord, hatMc, airVectorSub);
fDGbac      = projectAlg2DGsub(g, bacF, p, ord, hatMc, airVectorSub);
K11DG       = projectAlg2DGsub(g, K11, p, ord, hatMc, airVectorSub);
K12DG       = projectAlg2DGsub(g, K12, p, ord, hatMc, airVectorSub);
K21DG       = projectAlg2DGsub(g, K21, p, ord, hatMc, airVectorSub);
K22DG       = projectAlg2DGsub(g, K22, p, ord, hatMc, airVectorSub);
oneDG       = projectAlg2DGsub(g, one, p, ord, hatMc, airVectorSub);

NT = numAir;
N = (p+1)^2;

globL       = globM * reshape(fDG', NT*N, 1);
globLbac    = globM * reshape(fDGbac', NT*N, 1);

globG       = assembleGlobGsub(g, hatGc, hatGx, hatGy, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
globR       = assembleGlobRsub(g, markE0TairySub, hatRdiag, hatRoffdiag, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
globKN      = assembleGlobKNsub(g, p, ord, markE0TneumSub, g_N, airVectorSub, numAir);
globKNbac   = assembleGlobKNsub(g, p, ord, markE0TneumSub, bacG_N, airVectorSub, numAir);

sysOne      = reshape(oneDG', NT*N, 1);

sysA = [    globM                   ,   sparse(NT*N,NT*N)       ,   globH{1} + globQ{1} + globQN{1}     ;
            sparse(NT*N,NT*N)       ,   globM                   ,   globH{2} + globQ{2} + globQN{2}     ;
            globG{1} + globR{1}     ,   globG{2} + globR{2}     ,   globS                               ];
sysV1= [    -sysOne * g.deltaX^2    ;   zeros(size(globM,1),1)  ;   globKN + globL                      ];
sysV2= [    zeros(size(globM,1),1)  ;   -sysOne * g.deltaX^2    ;   globKN + globL                      ];
sysVbac = [ zeros(size(globM,1),1)  ;   zeros(size(globM,1),1)  ;   globKNbac + globLbac                ];
  
sysX1   = sysA \ sysV1;
sysX2   = sysA \ sysV2;
sysXbac = sysA \ sysVbac;

K11DG       = projectAlg2DGsub(g, K11, p, ord, hatMc, ones(g.numT,1));
K12DG       = projectAlg2DGsub(g, K12, p, ord, hatMc, ones(g.numT,1));
K21DG       = projectAlg2DGsub(g, K21, p, ord, hatMc, ones(g.numT,1));
K22DG       = projectAlg2DGsub(g, K22, p, ord, hatMc, ones(g.numT,1));

uDG                     = zeros(g.numTsub, N);
uDG(airVectorSub,:)     = reshape( sysX1( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
uLagr                   = projectDG2LagrangeSub(uDG);
visualizeDataSub(g, uLagr, 'u', 'u', 1);
uDG                     = zeros(g.numTsub, N);
uDG(airVectorSub,:)     = reshape( sysX2( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';
uLagr                   = projectDG2LagrangeSub(uDG);
visualizeDataSub(g, uLagr, 'u', 'u', 2);
q1DG                    = zeros(g.numTsub, N);
q1DG(airVectorSub,:)    = -reshape( sysX1( 0 * NT * N + 1 : 1 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta1_1                = reshape(q1DG', g.numT*N, 1);
q1Lagr                  = projectDG2LagrangeSub(q1DG);
visualizeDataSub(g, q1Lagr, 'theta1', 'theta1', 1);
q1DG                    = zeros(g.numTsub, N);
q1DG(airVectorSub,:)    = -reshape( sysX2( 0 * NT * N + 1 : 1 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta2_1                = reshape(q1DG', g.numT*N, 1);
q1Lagr                  = projectDG2LagrangeSub(q1DG);
visualizeDataSub(g, q1Lagr, 'theta1', 'theta1', 2);
q2DG                    = zeros(g.numTsub, N);
q2DG(airVectorSub,:)    = -reshape( sysX1( 1 * NT * N + 1 : 2 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta1_2                = reshape(q2DG', g.numT*N, 1);
q2Lagr                  = projectDG2LagrangeSub(q2DG);
visualizeDataSub(g, q2Lagr, 'theta2', 'theta2', 1);
q2DG                    = zeros(g.numTsub, N);
q2DG(airVectorSub,:)    = -reshape( sysX2( 1 * NT * N + 1 : 2 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta2_2                = reshape(q2DG', g.numT*N, 1);
q2Lagr                  = projectDG2LagrangeSub(q2DG);
visualizeDataSub(g, q2Lagr, 'theta2', 'theta2', 2);

visualizeDataSub(g, bulkVector, 'distribution', 'distribution', 1)

K_hom       = assembleGlobInt(g, hatInt, K11DG, K12DG, K21DG, K22DG, theta1_1, theta1_2, theta2_1, theta2_2)


uDGbac      = zeros(g.numTsub, N);
uDGbac(airVectorSub,:) = reshape( sysXbac( 2 * NT * N + 1 : 3 * NT * N ), N, NT )';

toc

end