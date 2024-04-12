function MainAggregat(  )

tic;

%% Definition of parameters                                                 %% Can be changed
NX          = 10;                       % Number of horicontal Elements
NZd         = 10;                       % Number of vertical Elements

tau         = 1;                        % End-Time
epsilon     = 1;                        % Everything < epsilon is defined as zero
numInnerIt  = 10;                       % Number of Diffusion/Absorption Iterations
numOuterIt  = 20;                       % Number of Transport/Growth Iterations

p           = 1;                        % Order of Polynomials for LDG Disc.
ord         = 4;                        % Order of Gaussian Integration Rule
eta         = 1;                        % Penalty Parameter for LDG

startConBio = 90;                       % Initial Concentration of Biomass
maxConcBio  = 100;                      % Max. Concentration of Biomass -> Growth
minConcBio  = 5;                        % Min. Concentration of Biomass -> Shrinkage

f_uptakeAgent = 20;                     % Growth Rate of Agent
f_uptakeAir = -epsilon/2;               % Amount of consumed Air when Agent grows
f_decay     = -10;                      % Degeneration Rate of Agent
f_bioMgrow  = 20;                       % Growth Rate of Biomass
f_bioMdecr  = -1;                       % Degeneration Rate of Biomass
f_bioUptAir = -epsilon/3;               % Amount of Air consumed, when Agent grows

%% Fixed Parameters (of general framework)
NZu         = 0;
width       = NX;
intBound    = NZd * ones(NX+1,1);
upBound     = NZd * ones(NX+1,1);

zero        = @(x,y) x-x + 0;
one         = @(x,y) x-x + 1;
uZero       = @(x,y) 100 / NX^2 * x .* (NX-x);
bacZero     = @(x,y) 100 / NX^2 * x .* (NX-x);
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
g = createDomain(width, NX, NZd, NZu, intBound, upBound);
markEbdr        = g.idE == -1 | g.idE == 1 | g.idE == 2 | g.idE == 6 | g.idE == 4;
concAgent       = 10*ones( g.numE , 1 );    % Initial Concentration of Agent
NX = g.NX;

%% Creating Initial Distribution of Bulk and Biomass                        %% Can be changed
% bulkVector      = randi([0,1], g.numTsub, 1); % 1 is bulk, 0 is air     % Vector containing Bulk Distribution

bulkVector = [zeros(43,1); ones(4,1); zeros(6,1); ones(4,1); zeros(6,1); ones(4,1); zeros(33,1)];
bioMvector = zeros(100,1);

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
bioConcVector   = startConBio * bioMvector;
biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', 0);
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
    
    %% Creating Vectors characterizing the Type of Edges
    idE             = zeros(g.numE, 1);
    idE(g.E0T(bulkVector == 1,1)) = idE(g.E0T(bulkVector == 1,1)) + 1;
    idE(g.E0T(bulkVector == 1,2)) = idE(g.E0T(bulkVector == 1,2)) + 1;
    idE(g.E0T(bulkVector == 1,3)) = idE(g.E0T(bulkVector == 1,3)) + 1;
    idE(g.E0T(bulkVector == 1,4)) = idE(g.E0T(bulkVector == 1,4)) + 1;
    idE = idE + markEbdr;

    idEneum         = idE == 1;
    idEairy         = idE == 0;
%     idEbulk         = ~(idEneum | idEairy);
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
        visualizeDataSub(g, uLagr, 'u', 'solu', 0);
        visualizeDataSub(g, bacLagr, 'bac', 'solBac', 0);
    end  % if

    %% Inner Loop for Diffusion and Growth/Shrinkage, ... (without Movement)
    
    for i = 1 : numInnerIt
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
        if i == 1 || i == numInnerIt
            uLagr       = projectDG2LagrangeSub( uDG );
            bacLagr     = projectDG2LagrangeSub( uDGbac );
            visualizeDataSub(g, uLagr, 'u', 'solu', (k-1)*numInnerIt+i);
            visualizeDataSub(g, bacLagr, 'bac', 'solBac', (k-1)*numInnerIt+i);
        end
        if i == 1
            biomVecLagr     = projectDG2LagrangeSub(bioConcVector);
            visualizeDataSub(g, biomVecLagr, 'bioMass', 'bioMass', (k-1)*numInnerIt+i);
        end
        
        concAir = computeConcAir(g, uDG, ord);
        concBac = computeConcAir(g, uDGbac, ord);
        
        %% Changes in Concentration of Agent        
        airyEdges = concAir(g.T0E(:,1)) > epsilon | concAir(g.T0E(:,2)) > epsilon;
        airyTrapezoids = zeros(g.numT,1);
        airyTrapezoids(g.T0E(airyEdges, 1)) = 1;
        airyTrapezoids(g.T0E(airyEdges, 2)) = 1;
        airyTrapezoids = airyTrapezoids .* (concAir > epsilon);
        airyTrapezoids = airyTrapezoids > 0;
        concAgent(airyEdges) = concAgent(airyEdges) + tau * (f_uptakeAgent + f_decay * concAgent(airyEdges));
        concAgent(~airyEdges) = concAgent(~airyEdges) + tau * f_decay * concAgent(~airyEdges);
        concAgent = concAgent .* (concAgent > 0);
        uDG(airyTrapezoids, 1) = uDG(airyTrapezoids, 1) + tau * f_uptakeAir;
        %% Changes in Concentration of Biomass
        bioMvector = bioMvector | (concBac > minConcBio);
        bioConcVector(bioMvector) = max(bioConcVector(bioMvector), minConcBio);
        bioMair = bioMvector .* (concAir > epsilon);
        bioMair = bioMair > 0;
        bioMnoAir = bioMvector .* (concAir <= epsilon);
        bioMnoAir = bioMnoAir > 0;
        bioConcVector(bioMair) = bioConcVector(bioMair) + tau * f_bioMgrow;
        bioConcVector(bioMnoAir) = bioConcVector(bioMnoAir) + tau * f_bioMdecr;
        bioMvector = bioConcVector > minConcBio;
        bioConcVector = bioConcVector .* bioMvector;
        uDG(bioMair, 1) = uDG(bioMair, 1) + tau * f_bioUptAir;
        
%        concAir = computeConcAir(g, uDG, ord);
%        concAir = concAir < epsilon;
%        uDG(concAir,:) = 0;
    end  % for i

    q1DG        = zeros(g.numTsub, N);
    q1DG(airVectorSub,:) = reshape( sysX( 1 : NT * N ), N, NT )';
    q2DG        = zeros(g.numTsub, N);
    q2DG(airVectorSub,:) = reshape( sysX( NT * N + 1 : 2 * NT * N ), N, NT )';
    q1DGbac        = zeros(g.numTsub, N);
    q1DGbac(airVectorSub,:) = reshape( sysXbac( 1 : NT * N ), N, NT )';
    q2DGbac        = zeros(g.numTsub, N);
    q2DGbac(airVectorSub,:) = reshape( sysXbac( NT * N + 1 : 2 * NT * N ), N, NT )';
    
    %% Doing the Movement of Biomass  
    bioOverHead = bioConcVector - maxConcBio;
    OverHead    = bioOverHead > 0;
    bioOverHead = bioOverHead .* OverHead;
    bioConcVector(OverHead) = maxConcBio;
    
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
    
    [~, indMax] = max(values);
    aims2 = zeros(size(indMax'));
    for i = 1 : length(indMax)
        aims2(i) = aims(indMax(i),i);
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
    bioConcVector(aims) = bioConcVector(aims) + bioOverHead(candidates);
    bioOverHead = bioConcVector - maxConcBio;
    OverHead    = bioOverHead > 0;
    bioConcVector(OverHead) = maxConcBio;
    bioMvector = bioConcVector > minConcBio;

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
            if aims(i,j) > 0
                aims(i,j) = aims(i,j) .* airVector(aims(i,j)) .* ~bioMvector(aims(i,j));
            end
            if aims(i,j) > 0
                values(i,j) = concAgent(g.E0T(aims(i,j),1)) + concAgent(g.E0T(aims(i,j),2)) ...
                    + concAgent(g.E0T(aims(i,j),3)) + concAgent(g.E0T(aims(i,j),4)) ...
                    + idEneum(g.E0T(aims(i,j),1)) + idEneum(g.E0T(aims(i,j),2)) ...
                    + idEneum(g.E0T(aims(i,j),3)) + idEneum(g.E0T(aims(i,j),4));
            end
        end
    end

    [~, indMax] = max(values);
    aims2 = zeros(size(indMax'));
    for i = 1 : length(indMax)
        aims2(i) = aims(indMax(i),i);
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
    
end  % for k

toc

end