function [ bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, ...
    concAgent, edgeChargeVector, particleList ,flag] = moveTripleBulk_ladung( bulkSize , stencilLayers, ...
    g , bulkVector, bulkTypeVector, particleTypeVector, uDG, q1DG, q2DG, uDGbac, q1DGbac, q2DGbac, ...
    concAgent, edgeChargeVector, agent, markEbdr, N, bioMvector, NZd , fileID, particleList, ...
    sumAgent,typeflag,stencil_type, attraction_type, particle_index)


numBulkOld = sum( bulkVector );
rotation = 0;

bulkVector_before = bulkVector;
bulkTypeVector_before = bulkTypeVector;
particleTypeVector_before = particleTypeVector;
particleList_before = particleList;
uDG_before = uDG;
q1DG_before = q1DG;
q2DG_before = q2DG;
uDGbac_before = uDGbac;
q1DGbac_before = q1DGbac;
q2DGbac_before = q2DGbac;
concAgent_before = concAgent;



if rotation == 1 
%% Rotation: finding candidates
% fprintf( fileID , 'Rotation \n' );
% fprintf( fileID , 'bulk size: %d \n' , bulkSize );

% Gibt benachbarte Kästchen an 
% jede reihe entspricht einem kästchen.  unten/oben/rechts/links
bulk0T                       = zeros( g.numT , 4 );
bulk0T( g.NX + 1 : end , 1 ) = bulkVector( 1 : end-g.NX );
bulk0T( 1 : end - g.NX , 2 ) = bulkVector( g.NX + 1 : end );
bulk0T( 1 : end - 1 , 3 )    = bulkVector( 2 : end );
bulk0T( 2 : end , 4 )        = bulkVector( 1 : end - 1 ); 

% gibt an, wie viele Kästchen sich an einer Kante befinden
idE                                 = zeros( g.numE , 1 );
idE( g.E0T( bulkVector == 1 , 1 ) ) = idE( g.E0T( bulkVector == 1 , 1 ) ) + 1;
idE( g.E0T( bulkVector == 1 , 2 ) ) = idE( g.E0T( bulkVector == 1 , 2 ) ) + 1;
idE( g.E0T( bulkVector == 1 , 3 ) ) = idE( g.E0T( bulkVector == 1 , 3 ) ) + 1;
idE( g.E0T( bulkVector == 1 , 4 ) ) = idE( g.E0T( bulkVector == 1 , 4 ) ) + 1;
idE                                 = idE + markEbdr;

idEairy  = idE == 0;

if typeflag == 2 || typeflag == 3 %only goethit and illite are allowed to rotate
    candidates = particleList;
else
    candidates = [];
end

rotationDirection = zeros(size(candidates,1),1);

if ~isempty( candidates )   
     
%% Rotation: finding aims

    if particle_index == 1
         height = 1;
         width = 17;
         position1 =  findParticleInds(1,17);
         position2 =  findParticleInds(17,1);
     elseif particle_index == 2
         height = 2;
         width = 34;        
         position1 =  findParticleInds(2,34);
         position2 =  findParticleInds(34,2);
     elseif particle_index == 3
         height = 2;
         width = 6;        
         position1 =  findParticleInds(2,6);
         position2 =  findParticleInds(6,2);     
     elseif particle_index == 4
         height = 2;
         width = 30;
         position1 =  findParticleInds(2,30);
         position2 =  findParticleInds(30,2);       
     elseif particle_index == 5
         height = 8;
         width = 100;
         position1 =  findParticleInds(8,100);
         position2 =  findParticleInds(100,8); 
     elseif particle_index == 6
         height = 1;
         width = 1;
         position1 =  findParticleInds(1,1);
         position2 =  findParticleInds(1,1);  
     end

    % holds Signifigance of each aim depending on sticky agent and neighbours
    aimValue = zeros( 1 , 6 );
    % holds the 6 possible positions the bulk could be in after rotation
    help = zeros( 6 , size( candidates , 2 ) );
    % go throug all the candidates
    for i = 1 : size( candidates , 1 )      
        aimValue( : ) = 0;
        help( : , : ) = 0;
        % The first position is the current position ( no rotation )
        help( 1 , : ) = candidates( i , : );

        if  typeflag == 2 %bulkSize == 5
            if bulkSize == 3
                load('goethit_3.mat')
            elseif bulkSize == 5
                load('goethit_5.mat')
            elseif bulkSize == 15
                load('goethit_15.mat')
            else
                assert('wrong goethit size')
            end
            candStencil = stencil( g.NX , NZd , candidates(i,1) , 20 );
            help = [candidates(i,:); candStencil(goethit_hor); candStencil(goethit_vert);candidates(i,:) ];
            
        % typeflag = 3 in our scenario 
            
        elseif typeflag == 3
            candStencil = stencil( g.NX , NZd , candidates(i,1) , ceil(size(candidates,2)/2) + 1 );
 
           % doesn't work for 8x100 particle yet
            stenHelp = stencil(g.NX,NZd,candidates(i,1),4);
            if  (ismember(stenHelp(19),candidates(i,:))  && particle_index ~=20)|| (ismember(stenHelp(9),candidates(i,:)) && particle_index == 7) ...
                    || (ismember(stenHelp(33),candidates(i,:)) && particle_index == 20) % horizontal
            %if ismember(candidates(i,1)-1, candidates(i,:)) && ismember(candidates(i,1)+1, candidates(i,:))
                if particle_index == 1
                    height = 1;
                    width = 17;
                elseif particle_index == 2
                    height = 2;
                    width = 34;        
                elseif particle_index == 3
                    height = 2;
                    width = 6;                
                elseif particle_index == 4
                    height = 2;
                    width = 30;      
                elseif particle_index == 5
                    height = 8;
                    width = 100;
                elseif particle_index == 6
                    height = 1;
                    width = 1;  
                end

                centerAfterRotation = [findCenterAfterRotation(candidates(i,1), g.NX, height, width, 1); findCenterAfterRotation(candidates(i,1), g.NX, height, width, 2);
                     findCenterAfterRotation(candidates(i,1), g.NX, height, width, 3); findCenterAfterRotation(candidates(i,1), g.NX, height, width, 4)];
                
                 candStencilHorizontalRotation = stencil( g.NX , NZd , centerAfterRotation , ceil(size(candidates,2)/2) + 1);

                help = [candidates(i,:); candStencil(position2); candStencilHorizontalRotation(1,position2); candStencilHorizontalRotation(2,position2);...
                    candStencilHorizontalRotation(3,position2); candStencilHorizontalRotation(4,position2)]; 
    
            else 
                if particle_index == 1
                    height = 17;
                    width = 1;
                elseif particle_index == 2
                    height = 34;
                    width = 2;        
                elseif particle_index == 3
                    height = 6;
                    width = 2;                
                elseif particle_index == 4
                    height = 30;
                    width = 2;      
                elseif particle_index == 5
                    height = 100;
                    width = 8;
                elseif particle_index == 6
                    height = 1;
                    width = 1;
                end  

                
                centerAfterRotation = [findCenterAfterRotation(candidates(i,1), g.NX, height, width, 1); findCenterAfterRotation(candidates(i,1), g.NX, height, width, 2);
                     findCenterAfterRotation(candidates(i,1), g.NX, height, width, 3); findCenterAfterRotation(candidates(i,1), g.NX, height, width, 4)];
                
                 candStencilVerticalRotation = stencil( g.NX , NZd , centerAfterRotation , ceil(size(candidates,2)/2) + 1);


                help = [candidates(i,:); candStencil(position1); candStencilVerticalRotation(1,position1); candStencilVerticalRotation(2,position1);...
                    candStencilVerticalRotation(3,position1); candStencilVerticalRotation(4,position1)]; 
                
            end      
        end
        % compute significance of possible aims
        % go through the 6 possible positions after rotation
        for j  = 1 : 6    %1 = 0 deg, 2 = 90 deg, 3 = 180 deg, 4 = 270 deg
            if j == 2 || j == 3 || j == 6
                bulk_trafo = get_particle_trafo(candidates( i , : ), g, NZd, 1, particle_index);
            elseif j == 4 || j == 5
                bulk_trafo = get_particle_trafo(candidates( i , : ), g, NZd, 2, particle_index);
            else
                bulk_trafo = 1:size(candidates,2);
            end            
            for k = 1 : size( candidates , 2 )
                % exclude possitions where there is bulk already
                help( j , 1 ) = help( j , 1 ) * ( ~bulkVector( help( j , k ) ) || any( candidates( i , : ) == help( j , k ) ) );
            end % for k
            
            % old error, always help ( j , 1 ) == candidates( i , 1 ) !!!
%           if help( j , 1 ) ~= candidates( i , 1 ) && help( j , 1 ) ~= 0
            if help( j , 1 ) ~= 0
            % 
                for k = 1 : size( help , 2 )
                    
                    % stencil = 5 is weird
%                     neighbours = stencil( g.NX , NZd , help( j , k ) , 5 );
                    neighbours = stencil( g.NX , NZd , help( j , k ) , 1 );
                    chargesCandidate = edgeChargeVector( g.CE0T( candidates( i , bulk_trafo(k) ) , : ) );            
                    %
                    for l = 1 : 4
                       if l == 1
                           edgeCandidate = 1;
                           edgeNeighbour = 2;
                            if j == 2 || j == 3 || j == 6
                                edgeCandidate = 3;
                            elseif j == 4 || j == 5
                                edgeCandidate = 4;
                            end
                           m = 2;    % unterer Nachbar der freien Sprungstelle 
                       elseif l == 2 
                           edgeCandidate = 2;
                           edgeNeighbour = 1;
                           if j == 2 || j == 3 || j == 6
                                edgeCandidate = 4;
                           elseif j == 4 || j == 5
                                edgeCandidate = 3;
                           end
                           m = 5;    % oberer 
                       elseif l == 3
                           edgeCandidate = 3;
                           edgeNeighbour = 4;
                            if j == 2 || j == 3 || j == 6
                                edgeCandidate = 2;
                            elseif j == 4 || j == 5
                                edgeCandidate = 1;
                            end
                           m = 4;    % rechter 
                       else 
                           edgeCandidate = 4;
                           edgeNeighbour = 3;
                            if j == 2 || j == 3 || j == 6
                                edgeCandidate = 1;
                            elseif j == 4 || j == 5
                                edgeCandidate = 2;
                            end
                           m = 3;    % linker 
                       end 
                       switch attraction_type
                           case 1 % with volume charges (old)
                                aimValue( j ) = aimValue( j ) + (( concAgent( g.E0T( help( j , k ) , l ) ) * bulkVector( neighbours( 1 , m ) ) + 0*bulkVector( neighbours( 1 , m ) ) - 1*bulkTypeVector(candidates( i , 1 ))*bulkTypeVector(neighbours(1,m)) ) ...
                                                * ( ~any( candidates( i , : ) == neighbours( 1 , m ) )  ) );
                           case 2 % without charges
                                aimValue( j ) = aimValue( j ) + (( concAgent( g.E0T( help( j , k ) , l ) ) * bulkVector( neighbours( 1 , m ) ) + 1*bulkVector( neighbours( 1 , m ) ) - 0*bulkTypeVector(candidates( i , 1 ))*bulkTypeVector(neighbours(1,m)) ) ...
                                                * ( ~any( candidates( i , : ) == neighbours( 1 , m ) )  ) );                               
                           case 3 % with edgeCharges      
                                aimValue( j ) = aimValue( j ) - (chargesCandidate(edgeCandidate)*edgeChargeVector(g.CE0T(neighbours(1,m),edgeNeighbour))) * ( ~any( candidates( i , : ) == neighbours( 1 , m ) ));

                       end % switch attraction_type                 
                    end % for l 
                end % for k #
            end % if
        end % for j
        
        % find the maximum of each column and its corresponding index
        [maximo , indMax] = max( aimValue( : ) );
        % find aims with maximal value for each candidate ( might be more than
        % one )
        indices = find( maximo == aimValue( : ) );
        % if there is more than one aim with the same maximal value choose one
        % randomly
        chooser = randi( length( indices ) );
        if length( indices ) == 1 && indMax ~= indices( chooser ) 
            error( 'IDIOT' )
        end
        aims( i , : ) = help( indices( chooser ) , : );
        rotationDirection(i) = indices( chooser );
        
        % check if current candidate is connected with sticky agent to
        % other particle. if connected, no rotation
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
    rotationDirection = rotationDirection( aims( : , 1 ) > 0 , : );
    aims       = aims( aims( : , 1 ) > 0 , : );
    

end 




% fprintf( fileID , 'number of candidates with a valid aim: %d \n' , size( candidates , 1 ) );
    
if ~isempty( candidates )
   
%% Rotation: conflict resolution        

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
        rotationDirection = rotationDirection( aims( : , 1 ) > 0 , : );
        aims       = aims( aims( : , 1 ) > 0 , : );
        
    end % for i
end % if ~isempty( candidates )

if ~isempty( candidates )
    
    if( length( unique( aims ) ) ~= ( size( aims , 1 ) * size( aims , 2 ) ) )
        aims
    end
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
    rotationDirection = rotationDirection(aims( : , 1 ) ~= 0 );
    aims = aims( aims( : , 1 ) ~= 0 , : );
    
    
    for i = 1:size( candidates , 1 )
        assert( length( bulkVector( candidates( i , : ) ) ) == length( find( bulkVector( candidates( i , : ) ) ) ) , 'problems with conflict resolution rotation')
    end
end % if ~ isempty( candidates )   

% fprintf( fileID , 'number of candidates after conflict resolution: %d \n' , size( candidates , 1 ) );

%% Rotation: switch places of candidates and coresponding aims and sticky agent 

if ~isempty( candidates )
%      [candidates,aims] 
%     sumAgent = sum(concAgent);
    for i = 1 : size( aims , 1 )

%       find candidates, that are NO aims and switch them with aims that are
%       NO candidates
        aimsInd = [];
        candInd = [];
        for j = 1 : size( aims , 2 ) 
            if ~any( candidates( i , j ) == aims( i , : ) )
                candInd = [ candInd , j ];
            end
            if ~any( aims( i , j ) == candidates( i , : ) )
                aimsInd = [ aimsInd , j ];
            end
        end
            
         
        % set bulkVector of candidates to 0 ( means no bulk )
        bulkVector( candidates( i , : ) ) = 0;
        bulkType_tmp = bulkTypeVector(candidates(i , :));
        bulkTypeVector( candidates( i , :) ) = 0;
        particleType_tmp = particleTypeVector(candidates(i , :));
        particleTypeVector( candidates( i , :) ) = 0;

        % set bulkVector of aims to 1 ( means bulk )
        bulkVector( aims( i , : ) ) = 1;
        bulkTypeVector( aims( i, : ) ) = bulkType_tmp;
        particleTypeVector( aims( i, : ) ) = particleType_tmp;
        [~,loc] = ismember(candidates(i,:),particleList,'rows');
        particleList(loc,:) = aims(i,:);
%%     % Check if sticky agent has to be moved
%         assert(sumAgent == sum(concAgent),'sumAgent_rotation_before')
%         if ~isequal(candidates(i,:),aims(i,:)) % rotate sticky agent
%             [candidates(i,:);aims(i,:)]
%             if typeflag == 2 % Stäbchen
%                 tmpAgent = concAgent(g.E0T(candidates(i,:),:));
%                 concAgent( g.E0T( candidates( i , : ) , : )) = 0;
%                 for bulk = 1 : size(candidates(i,:),2)                                
%                     concAgent( g.E0T( aims(i, bulk ),:)) = concAgent( g.E0T( aims(i, bulk ),:)) +  [tmpAgent(bulk,4) tmpAgent(bulk,3) tmpAgent(bulk,2) tmpAgent(bulk,1)]';              
%                 end
%             elseif typeflag == 3
%                 if bulkSize == 6
%                     bulk_trafo = get_illite_trafo(candidates(i,:), g, NZd);
%                     tmpAgent = concAgent(g.E0T(candidates(i,:),:));
%                     concAgent( g.E0T( candidates( i , : ) , : )) = 0;
%                     for bulk = 1 : size(candidates(i,:),2) 
%                         concAgent(g.E0T(aims(i,bulk),:)) = concAgent(g.E0T(aims(i,bulk),:)) + [tmpAgent(bulk_trafo(bulk),3) tmpAgent(bulk_trafo(bulk),4) tmpAgent(bulk_trafo(bulk),2) tmpAgent(bulk_trafo(bulk),1)]';
%                     end
%                 elseif bulkSize == 24
%                 elseif bulkSize == 96
%                 end
%             end
%         end
%         assert(sumAgent == sum(concAgent),'sumAgent_rotation_after')


%%      %Move the charges on the edges
    if ~isequal(candidates(i,:),aims(i,:))
        if  rotationDirection(i) == 2 || rotationDirection(i) == 3 || rotationDirection(i) == 6
            bulk_trafo = get_particle_trafo(candidates( i , : ), g, NZd, 1, particle_index);
        elseif rotationDirection(i) == 4 || rotationDirection(i) == 5
            bulk_trafo = get_particle_trafo(candidates( i , : ), g, NZd, 2, particle_index);
        else
            bulk_trafo = 1:size(candidates,2);
        end 
         
        tmpCharges = edgeChargeVector( g.CE0T( candidates( i , bulk_trafo ) , : ) );
        edgeChargeVector( g.CE0T( candidates( i , : ) , : ) ) = 0;
        
        
        if  rotationDirection(i) == 2 || rotationDirection(i) == 3 || rotationDirection(i) == 6
            % gegen Uhrzeigersinn
            for bulk = 1 : size(candidates(i,:),2) 
                edgeChargeVector(g.CE0T(aims(i,bulk),:)) = [tmpCharges(bulk_trafo(bulk),4) tmpCharges(bulk_trafo(bulk),3) tmpCharges(bulk_trafo(bulk),1) tmpCharges(bulk_trafo(bulk),2)];
            end 
        elseif rotationDirection(i) == 4 || rotationDirection(i) == 5
            % im Uhrzeigersinn
            for bulk = 1 : size(candidates(i,:),2) 
                edgeChargeVector(g.CE0T(aims(i,bulk),:)) = [tmpCharges(bulk_trafo(bulk),3) tmpCharges(bulk_trafo(bulk),4) tmpCharges(bulk_trafo(bulk),2) tmpCharges(bulk_trafo(bulk),1)];
            end 
        else
            for bulk = 1 : size(candidates(i,:),2) 
                edgeChargeVector(g.CE0T(aims(i,bulk),:)) = [tmpCharges(bulk_trafo(bulk),1) tmpCharges(bulk_trafo(bulk),2) tmpCharges(bulk_trafo(bulk),3) tmpCharges(bulk_trafo(bulk),4)];
            end 
        end 
    end       
        
%%
        assert( numBulkOld == sum(bulkVector) , 'solid Bloecke stimmen nicht' )
        
        uDG( candidates( i , candInd ) , : )     = uDG( aims( i , aimsInd ) , : );
        q1DG( candidates( i , candInd ) , : )    = q1DG( aims( i , aimsInd ) , : );
        q2DG( candidates( i , candInd ) , : )    = q2DG( aims( i , aimsInd ) , : );
        uDGbac( candidates( i , candInd ) , : )  = uDGbac( aims( i , aimsInd ) , : );
        q1DGbac( candidates( i , candInd ) , : ) = q1DGbac( aims( i , aimsInd ) , : );
        q2DGbac( candidates( i , candInd ) , : ) = q2DGbac( aims( i , aimsInd ) , : );

%       set aims to zero
        uDG( aims( i , : ) , : )     = zeros( size( aims , 2 ) , N );
        q1DG( aims( i , : ) , : )    = zeros( size( aims , 2 ) , N );
        q2DG( aims( i , : ) , : )    = zeros( size( aims , 2 ) , N );
        uDGbac( aims( i , : ) , : )  = zeros( size( aims , 2 ) , N );
        q1DGbac( aims( i , : ) , : ) = zeros( size( aims , 2 ) , N );
        q2DGbac( aims( i , : ) , : ) = zeros( size( aims , 2 ) , N );

    end % for i
end % if ~isempty( candidates )


% fprintf( fileID , 'number of candidates after conflict resolution: %d \n \n' , size( candidates , 1 ) );
end %rotation
assert(numBulkOld == sum(bulkVector), 'solid Bloecke stimmen nicht')


% display('change flag=0');
flag = 0;

time1 = tic;

movement = 1;
if movement == 1
%% Movement: finding candidates

% fprintf( fileID , 'Movement \n' );
% fprintf( fileID , 'bulk size: %d \n' , bulkSize );

bulk0T                        = zeros( g.numT , 4 );
bulk0T( g.NX + 1 : end , 1 )  = bulkVector( 1 : end-g.NX );
bulk0T( 1 : end - g.NX , 2 )  = bulkVector( g.NX + 1 : end );
bulk0T( 1 : end - 1 , 3 )     = bulkVector( 2 : end );
bulk0T( 2 : end , 4 )         = bulkVector( 1 : end - 1 );

%
idE                                 = zeros( g.numE , 1 );
idE( g.E0T( bulkVector == 1 , 1 ) ) = idE( g.E0T( bulkVector == 1 , 1 ) ) + 1;
idE( g.E0T( bulkVector == 1 , 2 ) ) = idE( g.E0T( bulkVector == 1 , 2 ) ) + 1;
idE( g.E0T( bulkVector == 1 , 3 ) ) = idE( g.E0T( bulkVector == 1 , 3 ) ) + 1;
idE( g.E0T( bulkVector == 1 , 4 ) ) = idE( g.E0T( bulkVector == 1 , 4 ) ) + 1;
idE                                 = idE + markEbdr;

idEairy  = idE == 0;

candidates = particleList;
% fprintf('Size of particle: %d ', length(candidates))
time2 = tic;
if ~isempty( candidates ) 
    
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
        %aim( : ) = 0;
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
            bulkTypeVector=bulkTypeVector_before;
            particleTypeVector=particleTypeVector_before;
            particleList=particleList_before;
            uDG=uDG_before;
            q1DG=q1DG_before;
            q2DG=q2DG_before;
            uDGbac=uDGbac_before;
            q1DGbac=q1DGbac_before;
            q2DGbac=q2DGbac_before;
            concAgent=concAgent_before;
            return 
            
        end
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

conflictResolution = 1;
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
%         assert(sumAgent == sum(concAgent),'sumAgent_movement_before')
% 
%         tmpAgent = concAgent(g.E0T(candidates(i,:),:));
%         concAgent(g.E0T(candidates(i,:),:)) = 0;
%         concAgent(g.E0T(aims(i,:),:)) = concAgent(g.E0T(aims(i,:),:)) + tmpAgent;
%         assert(sum(concAgent<0)==0,'negativer agent')
%         assert(sumAgent == sum(concAgent),'sumAgent_movement_after')
        
        % set bulkVector of candidates to 0 ( means no bulk )
        bulkVector( candidates( i , : ) ) = 0;
        bulkType_tmp = bulkTypeVector(candidates(i , :));
        bulkTypeVector( candidates( i , :) ) = 0;
        particleType_tmp = particleTypeVector(candidates(i , :));
        particleTypeVector( candidates( i , :) ) = 0;
        
        % set bulkVector of aims to 1 ( means bulk )
        bulkVector( aims( i , : ) ) = 1;
        bulkTypeVector( aims( i, : ) ) = bulkType_tmp;
        particleTypeVector( aims( i, : ) ) = particleType_tmp;
        [~,loc] = ismember(candidates(i,:),particleList,'rows');
        particleList(loc , : ) = aims(i , : );
        
        % Move the charges on the edges
        tmpCharges = edgeChargeVector( g.CE0T( candidates( i , : ) , : ) );
        edgeChargeVector( g.CE0T( candidates( i , : ) , : ) ) = 0;
        edgeChargeVector(g.CE0T(aims(i, :),:)) = tmpCharges;


        assert(numBulkOld == sum(bulkVector), 'solid Bloecke stimmen nicht, movement')

       
        % set values of candidates to values of aims
        uDG( candidates( i , candInd ) , : )     = uDG( aims( i , aimsInd ) , : );
        q1DG( candidates( i , candInd ) , : )    = q1DG( aims( i , aimsInd ) , : );
        q2DG( candidates( i , candInd ) , : )    = q2DG( aims( i , aimsInd ) , : );
        uDGbac( candidates( i , candInd ) , : )  = uDGbac( aims( i , aimsInd ) , : );
        q1DGbac( candidates( i , candInd ) , : ) = q1DGbac( aims( i , aimsInd ) , : );
        q2DGbac( candidates( i , candInd ) , : ) = q2DGbac( aims( i , aimsInd ) , : );

        % set values of aims to zero
        uDG( aims( i , : ) , : )     = zeros( size( aims , 2 ) , N );
        q1DG( aims( i , : ) , : )    = zeros( size( aims , 2 ) , N );
        q2DG( aims( i , : ) , : )    = zeros( size( aims , 2 ) , N );
        uDGbac( aims( i , : ) , : )  = zeros( size( aims , 2 ) , N );
        q1DGbac( aims( i , : ) , : ) = zeros( size( aims , 2 ) , N );
        q2DGbac( aims( i , : ) , : ) = zeros( size( aims , 2 ) , N );
        
    end % for i 
%     fprintf( fileID , 'sum( concAgent ) : %f \n \n' , sum( concAgent ) ); 
end % if ~isempty( candidates )
flag = 0;
% fprintf('Time for movement: %d \n', toc(time4))

% fprintf('Time for jumping of one aggregate: %d \n', toc(time1))
end
