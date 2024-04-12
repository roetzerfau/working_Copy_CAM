function K_hom = SimpleHomogenizationFunc(  bulkVector )

tic;

%% Definition of parameters                                                 %% Can be changed
NX          = sqrt(length(bulkVector));                       % Number of horicontal Elements
NZd         = NX;                       % Number of vertical Elements

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
K11         = one;
K12         = zero;
K21         = zero;
K22         = one;
g_N         = zero;

%% Creating domain for Simulation
g = createDomainFolded(width, NX, NZd, NZu, intBound, upBound);
markEbdr        = g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8 | g.idE == 8;

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
K11DG       = projectAlg2DGsub(g, K11, p, ord, hatMc, airVectorSub);
K12DG       = projectAlg2DGsub(g, K12, p, ord, hatMc, airVectorSub);
K21DG       = projectAlg2DGsub(g, K21, p, ord, hatMc, airVectorSub);
K22DG       = projectAlg2DGsub(g, K22, p, ord, hatMc, airVectorSub);
oneDG       = projectAlg2DGsub(g, one, p, ord, hatMc, airVectorSub);

NT = numAir;
N = (p+1)^2;

globL       = globM * reshape(fDG', NT*N, 1);

globG       = assembleGlobGsub(g, hatGc, hatGx, hatGy, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
globR       = assembleGlobRsub(g, markE0TairySub, hatRdiag, hatRoffdiag, K11DG, K12DG, K21DG, K22DG, airVectorSub, numAir);
globKN      = assembleGlobKNsub(g, p, ord, markE0TneumSub, g_N, airVectorSub, numAir);

sysOne      = reshape(oneDG', NT*N, 1);

sysA = [    globM                   ,   sparse(NT*N,NT*N)       ,   globH{1} + globQ{1} + globQN{1}     ;
            sparse(NT*N,NT*N)       ,   globM                   ,   globH{2} + globQ{2} + globQN{2}     ;
            globG{1} + globR{1}     ,   globG{2} + globR{2}     ,   globS                               ];
sysV1= [    -sysOne * g.deltaX^2    ;   zeros(size(globM,1),1)  ;   globKN + globL                      ];
sysV2= [    zeros(size(globM,1),1)  ;   -sysOne * g.deltaX^2    ;   globKN + globL                      ];
  
sysX1   = sysA \ sysV1;
sysX2   = sysA \ sysV2;

K11DG       = projectAlg2DGsub(g, K11, p, ord, hatMc, ones(g.numT,1));
K12DG       = projectAlg2DGsub(g, K12, p, ord, hatMc, ones(g.numT,1));
K21DG       = projectAlg2DGsub(g, K21, p, ord, hatMc, ones(g.numT,1));
K22DG       = projectAlg2DGsub(g, K22, p, ord, hatMc, ones(g.numT,1));

q1DG                    = zeros(g.numTsub, N);
q1DG(airVectorSub,:)    = -reshape( sysX1( 0 * NT * N + 1 : 1 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta1_1                = reshape(q1DG', g.numT*N, 1);
q1DG                    = zeros(g.numTsub, N);
q1DG(airVectorSub,:)    = -reshape( sysX2( 0 * NT * N + 1 : 1 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta2_1                = reshape(q1DG', g.numT*N, 1);
q2DG                    = zeros(g.numTsub, N);
q2DG(airVectorSub,:)    = -reshape( sysX1( 1 * NT * N + 1 : 2 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta1_2                = reshape(q2DG', g.numT*N, 1);
q2DG                    = zeros(g.numTsub, N);
q2DG(airVectorSub,:)    = -reshape( sysX2( 1 * NT * N + 1 : 2 * NT * N ), N, NT )'; % Beware of the minus for homogenization flux
theta2_2                = reshape(q2DG', g.numT*N, 1);

K_hom       = assembleGlobInt(g, hatInt, K11DG, K12DG, K21DG, K22DG, theta1_1, theta1_2, theta2_1, theta2_2);

toc

end