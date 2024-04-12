function [ bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector, POMageVector, ...
     concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, particleList ,flag] = ...
    moveParticles( bulkSize , stencilLayers , g , bulkVector, bulkTypeVector, particleTypeVector, POMVector, POMconcVector,...
    POMageVector, concAgent, concPOMAgent, POMagentAge, edgeChargeVector, reactiveSurfaceVector, NZd , fileID, ...
    particleList ,sumAgent,typeflag,stencil_type, attraction_type, particle_index)

numBulkOld = sum( bulkVector );

bulkVector_before = bulkVector;
bulkTypeVector_before = bulkTypeVector;
POMVector_before = POMVector;
POMconcVector_before = POMconcVector;
POMageVector_before = POMageVector;
particleTypeVector_before = particleTypeVector;
particleList_before = particleList;

concAgent_before = concAgent;
concPOMAgent_before = concPOMAgent;
POMagentAge_before = POMagentAge;

% display('change flag=0');
flag = 0;

time1 = tic;

movement = 1;
if movement == 1
%% Movement: finding candidates

candidates = particleList;
% fprintf('Size of particle: %d ', length(candidates))
time2 = tic;
if ~isempty( candidates ) 
    solidPOMreactiveEdgeIndicator = 0;
    solidPOMmemoryEdgeIndicator = 0;
    
%% Movement: finding aims
    helping = zeros( size( candidates , 1 ) , 2 * stencilLayers * ( stencilLayers + 1 ) + 1 , size( candidates , 2 ) );
    % AnzKand x AnzNachbarn(stencil) x Kandgröße
    for k = 1 : size( candidates, 2 )
        helping( : , : , k ) = stencil( g.NX, NZd, candidates( : , k ) , stencilLayers );
    end
    
    % aims is the chosen aim, only chosen at the end
    aims = zeros( size( candidates , 1 ) , bulkSize );
    
    for i = 1 : size( candidates , 1 ) % iteriere über kandidaten

        % holds the significance of the possible positions
        aim = zeros( 1 , size( helping , 2) );
        %aim( : ) = 0;help
        for j = 1 : size( helping , 2 ) % iteriere über nachbarn des kandidaten
            for k = 1 : size( candidates , 2 )
                % exclude possitions where there is bulk already except for the candidate itself
                helping( i , j , 1 ) = helping( i , j , 1 ) * ( ( ~bulkVector( helping( i , j , k ) ) ) || ( any( candidates( i , : ) == helping( i , j , k ) ) ) );
            end % for k
            %if helping( i , j , 1 ) ~= candidates( i , 1 ) && helping( i , j , 1 ) ~= 0
            if helping( i , j , 1 ) ~= 0 %potentiell freie Sprungstelle
                for k = 1 : size( candidates , 2 ) %iteriere über Blöcke des Kandidaten
                   neighbours = stencil( g.NX , NZd , helping( i , j , k ) , 1 ); % Nachbarn der potentiell freien Sprungstelle
                   for l = 1 : 4
                       if l == 1
                           edgeCandidate = 1;
                           edgeNeighbour = 2;
                           m = 2;    % unterer Nachbar der freien Sprungstelle 
                       elseif l == 2 
                           edgeCandidate = 2;
                           edgeNeighbour = 1;
                           m = 5;    % oberer 
                       elseif l == 3
                           edgeCandidate = 3;
                           edgeNeighbour = 4;
                           m = 4;    % rechter 
                       else 
                           edgeCandidate = 4;
                           edgeNeighbour = 3;
                           m = 3;    % linker 
                       end  
                       % compute the significance of a single bulk cell
                       % sum up all the significance values of one bulk
                       switch attraction_type
                           case 1 % with volume charges (old)
                                aim( j ) = aim( j ) + (( concAgent( g.E0T( helping( i , j , k ) , l ) ) * bulkVector( neighbours( 1 , m ) ) + 0*bulkVector( neighbours( 1 , m ) ) - 1*bulkTypeVector(candidates( i , k )) * bulkTypeVector(neighbours(1,m)) ) ...
                                              * ( ~any( candidates( i , : ) == neighbours( 1 , m ) )));
                           case 2 % without charges
                                aim( j ) = aim( j ) + (( concAgent( g.E0T( helping( i , j , k ) , l ) ) * bulkVector( neighbours( 1 , m ) ) + 1*bulkVector( neighbours( 1 , m ) ) - 0*bulkTypeVector(candidates( i , k )) * bulkTypeVector(neighbours(1,m)) ) ...
                                              * ( ~any( candidates( i , : ) == neighbours( 1 , m ) )));                           
                           case 3 % with edgeCharges      
                                aim( j ) = aim( j ) - (edgeChargeVector(g.CE0T(candidates( i , k ),edgeCandidate))*edgeChargeVector(g.CE0T(neighbours(1,m),edgeNeighbour))) * ( ~any( candidates( i , : ) == neighbours( 1 , m ) ));     
                           case 4 % with edgeCharges    
                                 aim( j ) = aim( j ) + (bulkVector(candidates( i , k ))* POMVector(neighbours( 1 , m )) * edgeChargeVector(g.CE0T(candidates( i , k ),edgeCandidate)) ...
                                     + bulkVector(neighbours( 1 , m ))* POMVector(candidates( i , k ))*edgeChargeVector(g.CE0T(neighbours(1,m),edgeNeighbour))) * ( ~any( candidates( i , : ) == neighbours( 1 , m ) ));
                           case 5 % with edgeCharges    
%                                  solidSolidAttr = (bulkVector(candidates( i , k ))* bulkVector(neighbours( 1 , m )) * reactiveSurfaceVector(g.CE0T(candidates( i , k ),edgeCandidate)) ...
%                                      + bulkVector(neighbours( 1 , m ))* bulkVector(candidates( i , k ))*reactiveSurfaceVector(g.CE0T(neighbours(1,m),edgeNeighbour))) * ( ~any( candidates( i , : ) == neighbours( 1 , m ) ));
                                 solidSolidAttr = (bulkVector(candidates( i , k ))* bulkVector(neighbours( 1 , m )) * reactiveSurfaceVector(g.CE0T(candidates( i , k ),edgeCandidate)) ...
                                     *reactiveSurfaceVector(g.CE0T(neighbours(1,m),edgeNeighbour))) * ( ~any( candidates( i , : ) == neighbours( 1 , m ) ));
                                 solidPOMreactiveSurfAttr = (bulkVector(candidates( i , k ))* POMVector(neighbours( 1 , m )) * reactiveSurfaceVector(g.CE0T(candidates( i , k ),edgeCandidate)) ...
                                     + bulkVector(neighbours( 1 , m ))* POMVector(candidates( i , k ))*reactiveSurfaceVector(g.CE0T(neighbours(1,m),edgeNeighbour))) * ( ~any( candidates( i , : ) == neighbours( 1 , m ) ));
                                 solidPOMmemoryAttr = (bulkVector(candidates( i , k ))* POMVector(neighbours( 1 , m )) * edgeChargeVector(g.CE0T(candidates( i , k ),edgeCandidate)) ...
                                     + bulkVector(neighbours( 1 , m ))* POMVector(candidates( i , k ))*edgeChargeVector(g.CE0T(neighbours(1,m),edgeNeighbour))) * ( ~any( candidates( i , : ) == neighbours( 1 , m ) ));
                                 aim( j ) = aim( j ) + max([solidSolidAttr 5*solidPOMreactiveSurfAttr 10*solidPOMmemoryAttr]);     
                                 if(solidPOMreactiveSurfAttr > 0 && j == 1)
                                     solidPOMreactiveEdgeIndicator = 1;
                                 end
                                 if(solidPOMmemoryAttr > 0 && j == 1)
                                     solidPOMmemoryEdgeIndicator = 1;
                                 end
    
                       end % switch attraction_type                       
                   end % for l
                end % for k 
            else
                aim(j) = -1;
            end % helping( i , j , 1 ) ~= candidates( i , 1 ) && helping( i , j , 1 ) ~= 0
        end % for j 
        
        % find the maximum of each column and its corresponding index
        [maximo , indMax] = max( aim( : ) );
        if maximo == 0 && stencil_type == 1 && aim(1)==maximo
            flag = 1;
            bulkVector =bulkVector_before;
            POMVector = POMVector_before;
            POMconcVector = POMconcVector_before;
            POMageVector = POMageVector_before;
            bulkTypeVector=bulkTypeVector_before;
            particleTypeVector=particleTypeVector_before;
            particleList=particleList_before;
            concAgent=concAgent_before;
            concPOMAgent = concPOMAgent_before;
            POMagentAge = POMagentAge_before;
            return 
            
        end
        
        % enable random breaking up
        % if current position is the most attractive
%         if aim(1)==maximo
            if solidPOMmemoryEdgeIndicator == 1
               randNum = randi(100,1);
               if randNum > 95
                   aim(1) = -1;
                   aim( aim >= 0 ) = 0; 
                   [maximo , indMax] = max( aim( : ) );
               end
            elseif solidPOMreactiveEdgeIndicator == 1
               randNum = randi(100,1);
               if randNum > 90
                   aim(1) = -1;
                   aim( aim >= 0 ) = 0; 
                   [maximo , indMax] = max( aim( : ) );
               end
            else
               randNum = randi(100,1);
               if randNum > 75
%                    aim(1) = -1;
                   aim( aim >= 0 ) = 0; 
                   [maximo , indMax] = max( aim( : ) );
               end
            end
        
        end
%         end
        
        % find aims with maximal value for each candidate ( might be more than one )
        indices = find( maximo == aim( : ) );


        % if there is more than one aim with the same maximal value choose one randomly
        chooser = randi( length( indices ) );
        if length( indices ) == 1 && indMax ~= indices( chooser ) 
            error( 'IDIOT' )
        end       
        aims( i , : ) = helping( i , indices( chooser ) , : );
        
        
%         possible_neighbours = stencil(g.NX,NZd,candidates(i,:),1);
%         possible_neighbours = reshape(possible_neighbours,1,size(possible_neighbours,1)*size(possible_neighbours,2));
%         possible_neighbours = setdiff(possible_neighbours,candidates(i,:),'legacy');
%         neighbours = possible_neighbours( bulkVector(possible_neighbours) == 1 );
%         neighbours_edges = unique( reshape( g.E0T(neighbours,:),1,size(g.E0T(neighbours,:),1)*size(g.E0T(neighbours,:),2) ) );
%         candidate_edges = unique( reshape( g.E0T(candidates(i,:),:),1,size(g.E0T(candidates(i,:),:),1)*size(g.E0T(candidates(i,:),:),2) ) );
%         common_edges = intersect(neighbours_edges,candidate_edges);
%         if sum(concAgent(common_edges))~=0 
%            aims(i,:)=candidates(i,:);
%         end
  
    end % for i 

    candidates = candidates( aims( : , 1 ) > 0 , : );
    aims       = aims( aims( : , 1 ) > 0 , : );
end % if ~isempty( candidates )

% fprintf( fileID , 'number of candidates with a valid aim: %d \n' , size( candidates , 1 ) );
% fprintf('Time for finding aims: %d \n', toc(time2))

conflictResolution = 0;
time3 = tic;
if ~isempty( candidates ) && conflictResolution

%% Movement: Conflict resolution

    % find cells that are the aim of multiple candidates 
    for i = 1 : size( aims , 2 )
        
        % set all but one aim, that hold at least one equal cell to zero
        for j = 1 : size( aims , 1 )
            if aims( j , i ) ~= 0
                [ conflictIdx1 , conflictIdx2 ] = find( aims( : , i : end ) == aims( j , i ) );
                %assert( isempty( setdiff( conflictIdx1 , index ) ) , 'problem' )
                randConflict = conflictIdx1( randi( length( conflictIdx1 ) ) );
                conflicts = setdiff( conflictIdx1 , randConflict );
                for k = 1 : length( conflicts )
                    for l = i : size( aims , 2 )
                        aims( conflicts( k ) , : ) = aims( conflicts( k ) , : ) * ~any( aims( conflicts( k ) , : ) == aims( randConflict , l ) );
                    end % for l
                end % for k
            end % if
        end % for j
        
        % erase the conflicting aims and candidates
        candidates = candidates( aims( : , 1 ) > 0 , : );
        aims       = aims( aims( : , 1 ) > 0 , : );
    end % for i
end % if ~isempty( candidates )

if ~isempty( candidates ) && conflictResolution
    
    assert( length( unique( aims ) ) == ( size( aims , 1 ) * size( aims , 2 ) ) , 'there are still conflicts' )
    
    assert( length( candidates( : , 1 ) ) == length( aims( : , 1 ) ) , 'Laenge passt nicht' );
    
    neighbours = zeros( size( candidates , 1 ) , 5 , size( candidates , 2 ) );
    
    conflicts = zeros( size( candidates , 1 ) , 4 , size( candidates , 2 ) );
    
    for k = 1 : size( candidates , 2 )
        neighbours( : , : , k ) = stencil( g.NX , NZd , aims( : , k ) , 1 );
    end % for k
    neighbours = neighbours( : , 2 : 5 , : );
    
    for k = 1 : size( candidates , 2 ) 
        for n = 1 : 4 
            for i = 1 : length( aims( : , 1 ) ) 
            % a conflict occurs if the neighbour of an aim IS a bulk AND
            % the neighbour of an aim IS NOT the candidate of the aim AND
            % if the neighbour IS any candidate
                if bulkVector( neighbours( i , n , k ) ) == 1 && candidates( i , k ) ~= neighbours( i , n , k ) && any( candidates( : , k ) == neighbours( i , n , k ) )
                    conflicts( i , n , k ) = 1;
                end % if
            end % for i
        end % for n
    end % for k
    
    while sum( sum( sum( conflicts ) ) ) ~= 0
        % get all the neighbours that have conflicts
        problems                           = neighbours( conflicts == 1 );
        % get one random index
        index                              = randi( length( problems ) );
        % get the problematic neighbours at the random index
        mf                                 = problems( index );
        % find the candidate that is this problematic neighbour
        [ cand1 , cand2 ]                  = find( candidates( : , : ) == mf  );
        % set the problematic neigbour, its aims and its conflicts to 0
        aims( cand1 , : )                  = 0;
        candidates( cand1 , : )            = 0;
        conflicts( cand1 , : , : )         = 0;
        % 
        for j = 1 : 4
            [ cand1 , cand2 ] = find( neighbours( : , j , : ) == mf );
            conflicts( cand1 , j , cand2 ) = 0;
        end % for i 
    end % while

    candidates = candidates( aims( : , 1 ) ~= 0 , : );
    aims = aims( aims( : , 1 ) ~= 0 , : );
end % if ~isempty( candidates ) 
    
% fprintf( fileID , 'number of candidates after conflict resolution: %d \n \n' , size( candidates , 1 ) );
% fprintf('Time for conflict resolution: %d \n', toc(time3))
time4 = tic;
%% Movement: switch places of candidates and coresponding aims
% [candidates,aims]
if ~isempty( candidates )
    % move the bulk, candidates and aims are swaped
    
    for i = 1 : size( aims , 1 )
        
        % find candidates, that are NO aims and switch them with aims that are
        % NO candidates
        aimsInd = [];
        candInd = [];
        for j = 1 : size( aims , 2 ) 
            if ~any( candidates( i , j ) == aims( i , : ) )
                candInd = [ candInd , j ];
            end % if
            if ~any( aims( i , j ) == candidates( i , : ) )
                aimsInd = [ aimsInd , j ];
            end % if
        end % for j  
        assert( sum(bulkVector(candidates(i,:)))== bulkSize,'candidate not complete')

        % set bulkVector of candidates to 0 ( means no bulk )
        bulkVector( candidates( i , : ) ) = 0;
        bulkType_tmp = bulkTypeVector(candidates(i , :));
        bulkTypeVector( candidates( i , :) ) = 0;
        POM_tmp = POMVector(candidates(i , :));
        POMVector(candidates(i , :)) = 0;
        POMconc_temp = POMconcVector(candidates(i , :));
        POMconcVector(candidates(i , :)) = 0;
        POMageVector_temp = POMageVector(candidates(i , :));
        POMageVector(candidates(i , :)) = 0;
        particleType_tmp = particleTypeVector(candidates(i , :));
        particleTypeVector( candidates( i , :) ) = 0;
        
        % set bulkVector of aims to 1 ( means bulk )
        bulkVector( aims( i , : ) ) = 1;
        bulkTypeVector( aims( i, : ) ) = bulkType_tmp;
        POMVector( aims( i, : ) ) =  POM_tmp;
        POMconcVector( aims( i, : ) ) =  POMconc_temp;
        POMageVector( aims( i, : ) ) =  POMageVector_temp;
        particleTypeVector( aims( i, : ) ) = particleType_tmp;
        [~,loc] = ismember(candidates(i,:),particleList,'rows');
        particleList(loc , : ) = aims(i , : );
        
        % Move the charges and the agent on the edges
        tmpPOMAgent = concPOMAgent(g.CE0T( candidates( i , : ) , : ));
        concPOMAgent(g.CE0T( candidates( i , : ) , : )) = 0;
        concPOMAgent(g.CE0T(aims(i, :),:)) = tmpPOMAgent;
        
        tmpPOMagentAge = POMagentAge(g.CE0T( candidates( i , : ) , : ));
        POMagentAge(g.CE0T( candidates( i , : ) , : )) = 0;
        POMagentAge(g.CE0T(aims(i, :),:)) = tmpPOMagentAge;
                
        tmpCharges = edgeChargeVector( g.CE0T( candidates( i , : ) , : ) );
        edgeChargeVector( g.CE0T( candidates( i , : ) , : ) ) = 0;
        edgeChargeVector(g.CE0T(aims(i, :),:)) = tmpCharges;
        
        tmpReactiveSurface = reactiveSurfaceVector( g.CE0T( candidates( i , : ) , : ) );
        reactiveSurfaceVector( g.CE0T( candidates( i , : ) , : ) ) = 0;
        reactiveSurfaceVector(g.CE0T(aims(i, :),:)) = tmpReactiveSurface;


        assert(numBulkOld == sum(bulkVector), 'solid Bloecke stimmen nicht, movement')
        
    end % for i 
%     fprintf( fileID , 'sum( concAgent ) : %f \n \n' , sum( concAgent ) ); 
end % if ~isempty( candidates )
flag = 0;
% fprintf('Time for movement: %d \n', toc(time4))

% fprintf('Time for jumping of one aggregate: %d \n', toc(time1))

end
