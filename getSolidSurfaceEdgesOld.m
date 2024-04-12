function [particleSurfaceEdgeList] = getSolidSurfaceEdgesOld(g, particleList, particleTypeVector)

% fileID = fopen('Input/testSurface.txt','r');
% formatSpec = '%f';
% particleTypes = fscanf(fileID,formatSpec);
% fclose(fileID);
% 
% numSquares = size(particleTypes,1);
% sizeStruc = sqrt(size(particleTypes,1));

% g = createDomainFolded(sizeStruc,sizeStruc,sizeStruc,0,sizeStruc*ones(sizeStruc+1,1),sizeStruc*ones(sizeStruc+1,1));


% edgeChargeVector = 0*ones( g.numCE , 1);   
% particleTypeVector = particleTypes;



% particleList = createSolidParticleList(particleTypeVector);


particleSurfaceEdgeList = cell(length(particleList),1);

for i = 1 : length(particleList)
    edgeCandidates = [];
    j = 1;
    while (isempty(edgeCandidates))
        solidInd = particleList{i}(j);
        possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
        possibleNeighbors = possibleNeighbors(2:end);
        for neigh = 1 : 4
            if (neigh == 1)
                edgeDirection = 1;
                if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                    edgeCandidates = g.CE0T(solidInd, edgeDirection);
                end
            elseif (neigh == 2)
                edgeDirection = 4;
                if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                    edgeCandidates = g.CE0T(solidInd, edgeDirection);
                end
            elseif (neigh == 3)
                edgeDirection = 3;
                if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                    edgeCandidates = g.CE0T(solidInd, edgeDirection);
                end
            elseif (neigh == 4)
                edgeDirection = 2;
                if (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) )
                    edgeCandidates = g.CE0T(solidInd, edgeDirection);
                end
            end
        end
        j = j + 1;
    end
    vertices = g.V0CE(edgeCandidates, :);
    visitedVertices = vertices(1);
    visitedEdges = edgeCandidates;
    currVertex = vertices(2);
    particleSurfaceEdgeList{i} = [particleSurfaceEdgeList{i} edgeCandidates(1)];
    while(~isempty(edgeCandidates) && ~isempty(currVertex))
        edgeCandidates = [];
        [possibleEdges,~] = find(g.V0CE == currVertex);
        possibleEdges = possibleEdges(~ismember(possibleEdges,visitedEdges));
        for indEdge = 1 : length(possibleEdges)
           solidInd = mod(possibleEdges(indEdge), g.NX * g.NX);
           possibleNeighbors = stencil( g.NX, g.NX, solidInd, 1); 
           possibleNeighbors = possibleNeighbors(2:end);
           [~,edgeDirection] = find(g.CE0T(solidInd,:)==possibleEdges(indEdge));
            for neigh = 1 : 4
                if ( (neigh == 1) && (edgeDirection == 1) )
                    if ((particleTypeVector(solidInd) == i) && (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) ))
                        edgeCandidates = g.CE0T(solidInd, edgeDirection);
                    end
                elseif ( (neigh == 2) && (edgeDirection == 4) )
                   if ((particleTypeVector(solidInd) == i) && (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) ))
                        edgeCandidates = g.CE0T(solidInd, edgeDirection);
                    end
                elseif ( (neigh == 3) && (edgeDirection == 3) )
                    if ((particleTypeVector(solidInd) == i) && (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) ))
                        edgeCandidates = g.CE0T(solidInd, edgeDirection);
                    end
                elseif ( (neigh == 4) && (edgeDirection == 2) )
                    if ((particleTypeVector(solidInd) == i) && (particleTypeVector(possibleNeighbors(neigh)) ~= particleTypeVector(solidInd) ))
                        edgeCandidates = g.CE0T(solidInd, edgeDirection);
                    end
                end
            end
        end
        if(~isempty(edgeCandidates))
            particleSurfaceEdgeList{i} = [particleSurfaceEdgeList{i} edgeCandidates(1)];
            visitedVertices = [visitedVertices currVertex];
            vertices = g.V0CE(edgeCandidates, :);
            vertices = vertices(~ismember(vertices,visitedVertices));
            visitedEdges = [visitedEdges edgeCandidates];
            currVertex = vertices;
        end        
    end
    
end

% for i = 1 : length(particleSurfaceEdgeList)
%    edgeChargeVector(particleSurfaceEdgeList{i}) = 1;
% end

end
