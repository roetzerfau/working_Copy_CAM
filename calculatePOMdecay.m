 function [bulkVector, bulkTypeVector, POMVector, POMconcVector, POMageVector, concPOMAgent, edgeChargeVector, POMParticleList,...
     POMsolidEdgeList, sumExcessPOM, POMagentInput] = calculatePOMdecay(g, parameters, bulkVector, bulkTypeVector, POMVector, POMconcVector,...
     POMageVector, concPOMAgent, edgeChargeVector, reactiveSurfaceVector, POMParticleList, sumExcessPOM, POMagentInput)
    
    
    numFluidNeighVector = calculateNumFluidNeighbors(g, bulkVector, edgeChargeVector, POMVector, 2);
    
%     POMconcVector(POMVector == 1) = POMconcVector(POMVector == 1) - parameters.POMdecayRate;
    
    % here concPOMAgent is considered as an edgeChargeVector
    POMsolidEdgeList = calculatePOMsolidEdgeList(g, bulkVector, POMVector, POMParticleList);
    
    % update POM conc
    POMconcVectorOld = POMconcVector;
    
%     for i = 1 : length(POMParticleList)
%         if sum(edgeChargeVector(POMsolidEdgeList{i})) > 0
%             POMconcVector(POMParticleList{i}) = POMconcVector(POMParticleList{i}) - parameters.POMdecayRate * numFluidNeighVector(POMParticleList{i});
%         end
%     end

%     % only POM particles that are attached to reactive solid surface arePOMconcVector(POMParticleList{i})
%     % decaying
%     for i = 1 : length(POMParticleList)
%         if sum(reactiveSurfaceVector(POMsolidEdgeList{i})) > 0
%             POMconcVector(POMParticleList{i}) = POMconcVector(POMParticleList{i}) - parameters.POMdecayRate * numFluidNeighVector(POMParticleList{i});
%         end
%     end
    tau = 1;
    % only POM particles that are attached to reactive solid surface are
    % decaying
    for i = 1 : length(POMParticleList)
        if (sum(reactiveSurfaceVector(POMsolidEdgeList{i})) + sum(edgeChargeVector(POMsolidEdgeList{i})) > 0)...
                && (sum(numFluidNeighVector(POMParticleList{i})) > 0)
            concOld = sum(POMconcVector(POMParticleList{i}));
            particleDecayRate = parameters.POMdecayRate * sum(numFluidNeighVector(POMParticleList{i})) / ...
                (sum(numFluidNeighVector(POMParticleList{i})) + length(POMsolidEdgeList{i}));
            concNew = concOld * exp(- particleDecayRate * tau);
            concDiff = concOld - concNew;
            POMconcVector(POMParticleList{i}) = POMconcVector(POMParticleList{i}) - concDiff * ...
                numFluidNeighVector(POMParticleList{i}) / sum(numFluidNeighVector(POMParticleList{i}));
        end
    end
    
%     POMconcVector(POMVector == 1) = POMconcVector(POMVector == 1) - parameters.POMdecayRate * numFluidNeighVector(POMVector == 1);
    concPOMAgent(concPOMAgent > 0) = concPOMAgent(concPOMAgent > 0) - parameters.POMagentDecayRate * concPOMAgent(concPOMAgent > 0);
    
    % set POM conc. below threshold to zero
    POMconcVector(POMconcVector < parameters.POMminConc) = 0;
    
    POMdecayVector = POMconcVectorOld - POMconcVector;
    
    concPOMAgent(concPOMAgent < 0) = 0;
    
    % save amount of decayed POM
    decayedPOMfromParticle = zeros(length(POMParticleList),1);
    for i = 1 : length(POMParticleList)
        decayedPOMfromParticle(i) = sum(POMdecayVector(POMParticleList{i}));
    end
    
    POMagentInput = POMagentInput + sum(decayedPOMfromParticle) * parameters.carbonUseEfficiency;
    
%     % get all solid-POM edges of POM particles where decay happened
%     POMsolidEdgeList = cell(1,length(POMParticleList));
%     for i = 1 : length(POMParticleList)
%         edgeCandidates = g.E0T(POMParticleList{i},:);
%         edgeCandidates = reshape(edgeCandidates, [size(edgeCandidates,1)*4,1]);
%         for j = 1 : length(edgeCandidates)
%             if( (bulkVector(g.T0E(edgeCandidates(j),1)) == 1 &&  POMVector(g.T0E(edgeCandidates(j),1)) == 0 ...
%                     && POMVector(g.T0E(edgeCandidates(j),2)) == 1 ) || ( bulkVector(g.T0E(edgeCandidates(j),2)) == 1 ...
%                     && POMVector(g.T0E(edgeCandidates(j),2)) == 0 && POMVector(g.T0E(edgeCandidates(j),1)) == 1 ) )
%                 
%                 POMsolidEdgeList{i} = [POMsolidEdgeList{i}; edgeCandidates(j)];
%                 
%             end
%         end
%     end

    % update POMagent on solid-POM edges
    for i = 1 : length(POMParticleList)
        if ~isempty(POMsolidEdgeList{i})
        concPOMAgent(POMsolidEdgeList{i}) = concPOMAgent(POMsolidEdgeList{i}) + parameters.carbonUseEfficiency ...
            * decayedPOMfromParticle(i) / length(POMsolidEdgeList{i});    
        end
    end
    
    % spreading
    for i = 1 : length(POMsolidEdgeList)
        for edge = 1 : length(POMsolidEdgeList{i})
           if concPOMAgent(POMsolidEdgeList{i}(edge)) > parameters.POMagentMax
               excessPOM = concPOMAgent(POMsolidEdgeList{i}(edge)) - parameters.POMagentMax;
               concPOMAgent(POMsolidEdgeList{i}(edge)) = parameters.POMagentMax;
               visitedVertices = g.V0CE(POMsolidEdgeList{i}(edge),:);
               visitedEdges = POMsolidEdgeList{i}(edge);
               [edgeHelper, ~] = find(ismember(g.V0CE, visitedVertices));
               edgeHelper = edgeHelper(~ismember(edgeHelper,visitedEdges));
               edgeCandidates = [];
               for indEdge = 1 : length(edgeHelper)
                  solidInd = mod(edgeHelper(indEdge), g.NX * g.NX);
                  if ((bulkVector(solidInd) == 1) &&  (POMVector(solidInd) == 0))
                      possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
                      possibleNeighbors = possibleNeighbors(2:end); 
                      [~,edgeDirection] = find(g.CE0T(solidInd,:)==edgeHelper(indEdge));
                      for neigh = 1 : 4
                        if ( (neigh == 1) && (edgeDirection == 1) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        elseif ( (neigh == 2) && (edgeDirection == 4) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        elseif ( (neigh == 3) && (edgeDirection == 3) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        elseif ( (neigh == 4) && (edgeDirection == 2) )
                            if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                            end
                        end
                      end
                  end
               end
             
               while ( (excessPOM > 0) &&  ~isempty(edgeCandidates) )
                   % get capacity
                   aim = edgeCandidates(1);
                   edgeCandidates = edgeCandidates(2:end);
                   visitedEdges = [visitedEdges aim];
                   
                   aimCapacity = max(parameters.POMagentMax - concPOMAgent(aim), 0);
                   
                   % edgeCandidate has enough capacity to take up the
                   % excessPOM
                   if(aimCapacity >= excessPOM)
                       concPOMAgent(aim) = concPOMAgent(aim) + excessPOM;
                       excessPOM = 0;
                   % aim can not take up all the excess POM    
                   else
                       if aimCapacity > 0
                          concPOMAgent(aim) = parameters.POMagentMax;
                       end
                       excessPOM = excessPOM - aimCapacity;
                       
                       verticesCandidates = g.V0CE(aim, :);
                       verticesCandidates = verticesCandidates(~ismember(verticesCandidates,visitedVertices));
                       % verticesCandidates should never be more than one,
                       % because every vertex comes from one that has
                       % already been visited
                       if ~isempty(verticesCandidates)
                           [edgeHelper, ~] = find(ismember(g.V0CE, verticesCandidates));
                           edgeHelper = edgeHelper(~ismember(edgeHelper,visitedEdges));
                           
                           for indEdge = 1 : length(edgeHelper)
                              solidInd = mod(edgeHelper(indEdge), g.NX * g.NX);
                              if ((bulkVector(solidInd) == 1) &&  (POMVector(solidInd) == 0))
                                  possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
                                  possibleNeighbors = possibleNeighbors(2:end); 
                                  [~,edgeDirection] = find(g.CE0T(solidInd,:)==edgeHelper(indEdge));
                                  for neigh = 1 : 4
                                    if ( (neigh == 1) && (edgeDirection == 1) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    elseif ( (neigh == 2) && (edgeDirection == 4) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    elseif ( (neigh == 3) && (edgeDirection == 3) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    elseif ( (neigh == 4) && (edgeDirection == 2) )
                                        if ((POMVector(possibleNeighbors(neigh)) == 1) || (bulkVector(possibleNeighbors(neigh)) == 0))
                                            edgeCandidates = [edgeCandidates edgeHelper(indEdge)];
                                        end
                                    end
                                  end
                              end
                           end
                       end
                   end
                                   
                   
               end
               if (excessPOM > 0)
                    fprintf('POM agent could not be spread!')
                    sumExcessPOM = sumExcessPOM + excessPOM;
                    excessPOM = 0;
               end
%                assert(excessPOM == 0, 'POM agent could not be spread!')
               
           end
        end
    end
    
%     for i = 1 : length(POMsolidEdgeList)
%         indHelper = find(concPOMAgent(POMsolidEdgeList{i}) > 1);
%         
%         concPOMAgentHelper =  concPOMAgent(POMsolidEdgeList{i});
%         concPOMAgent(POMsolidEdgeList{i}(indHelper)) = 1;
%         
%         concPOMAgentDiff = concPOMAgentHelper - concPOMAgent(POMsolidEdgeList{i});
%         excessPOMAgent = sum(concPOMAgentDiff);
%         if(excessPOMAgent > 0)
%            % spread excessPOMAgent to neighboring edges  
%         end
%     end
    
    
%     % Update edgeChargeVector, find all candidates
%     newSolidCandidates = cell(1,length(POMParticleList));
%     newPOMNeighbors = cell(1,length(POMParticleList));
%     for i = 1 : length(POMParticleList)
% %         indHelper = find(~POMVector(g.T0E( POMsolidEdgeList{i},:)));
%         indHelper = POMVector(g.T0E( POMsolidEdgeList{i},:)) == 0;
% %         indHelperPOM = find(POMVector(g.T0E( POMsolidEdgeList{i},:)));
%         indHelperPOM = POMVector(g.T0E( POMsolidEdgeList{i},:)) == 1;
%         solidIndHelper =g.T0E( POMsolidEdgeList{i},:);
%         
%         solidIndHelper = solidIndHelper';
%         indHelper = indHelper';
%         indHelperPOM = indHelperPOM';
%         
%         newSolidCandidates{i} = solidIndHelper(indHelper);
%         newPOMNeighbors{i} = solidIndHelper(indHelperPOM);
%     end
%     
%     % find all charged edges where POM agent is positiv
%     for i = 1 : length(POMParticleList)
%         for j = 1 : length(newSolidCandidates{i})
%             if  concPOMAgent(POMsolidEdgeList{i}(j)) > 0
%                 
%                possibleNeighbors = stencil( g.NX, g.NX, newSolidCandidates{i}(j), 1); 
%                possibleNeighbors = possibleNeighbors(2:end); 
%                for neigh = 1 : 4
%                   if (possibleNeighbors(neigh) == newPOMNeighbors{i}(j))
%                      if (neigh == 1)
%                          edgeChargeVector(g.CE0T(newSolidCandidates{i}(j),1)) = 1;
%                      elseif (neigh == 2)
%                          edgeChargeVector(g.CE0T(newSolidCandidates{i}(j),4)) = 1;
%                      elseif (neigh == 3)
%                          edgeChargeVector(g.CE0T(newSolidCandidates{i}(j),3)) = 1;
%                      elseif (neigh == 4)
%                          edgeChargeVector(g.CE0T(newSolidCandidates{i}(j),2)) = 1;
%                      end
% 
%                   end
%                end 
%            end
%         end
%     end

%     concPOMAgent(concPOMAgent < parameters.POMagentMin) = 0;
    
    edgeChargeVector( concPOMAgent > parameters.POMagentMin ) = 1;


    % update bulkVector etc.
    indPOMchanged = find( ( POMVector == 1 ) &  ( POMconcVector == 0 ) );
    
    POMVector(indPOMchanged) = 0;
    bulkVector(indPOMchanged) = 0;
    bulkTypeVector(indPOMchanged) = 0;
    POMageVector(indPOMchanged) = 0;
    
    % Update POM age
    POMageVector(POMageVector > 0) = POMageVector(POMageVector > 0) + 1;

    
    % if POM cells dis                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             appeared
    if (~isempty(indPOMchanged))
%         [POMParticleList, ~] = createPOMParticleList(POMVector);
        % remove degraded POM cells from POMParticleList
        for i = 1 : length(POMParticleList)
            POMParticleList{i}(find(ismember(POMParticleList{i}, indPOMchanged))) = [];
        end
        % remove all POM particles that are completely degraded from POMParticleList
        for i = length(POMParticleList) : -1 : 1
            if (length(POMParticleList{i})) == 0
                POMParticleList(i) = [];
            end
        end
        
        % find POM particles that are separated into multiple parts due to
        % degradation
        POMParticleListOld = POMParticleList;
        POMparticlesToRemove = [];
        for i = 1 : length(POMParticleListOld)
            helperVector = zeros(size(POMVector));
            helperVector(POMParticleListOld{i}) = 1;
            % if particle was separated, length(POMListHelper) > 1
            % check periodicity...
            [POMListHelper, ~] = createPOMParticleListNew(helperVector);
            if length(POMListHelper) > 1
                % add separated POM particles at end and remove original
                % from POMParticleList
                POMParticleList = [POMParticleList POMListHelper];
                POMparticlesToRemove = [POMparticlesToRemove i];
            end
        end
        POMParticleList(POMparticlesToRemove) = [];
        
    end
    
    
 
 end
