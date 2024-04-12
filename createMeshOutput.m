function createMeshOutput( bulkVector , num )

fileName    = ['mesh.', num2str(num), '.geo'];
file        = fopen(fileName, 'wt');

N_all       = length(bulkVector);
NX          = sqrt(N_all);
h           = 1/NX;
N_air       = N_all - sum(bulkVector);

% numNoEdges  = 2 * sum( bulkVector(1:NX) + bulkVector(end-NX+1:end) == 2 ) ...
%                 + 2 * sum( bulkVector(1:NX:end) + bulkVector(NX:NX:end) == 2 );
%             
% for i = 0 : NX-1
%     numNoEdges = numNoEdges + sum( bulkVector(i*NX+1:(i+1)*NX-1) + bulkVector(i*NX+2:(i+1)*NX) == 2 ) ...
%                     + sum( bulkVector(i*NX+1:NX:end-NX) + bulkVector((i+1)*NX+1:NX:end) == 2 );
% end

N_edges = 0;

for i = 1 : NX
    N_edges = N_edges + ~bulkVector(i) + ~bulkVector((NX-1)*NX+i) + ~bulkVector((i-1)*NX+1) + ~bulkVector(i*NX) ...
        + sum( bulkVector((i-1)*NX+1:i*NX-1) + bulkVector((i-1)*NX+2:i*NX) == 1 ) ...
        + sum( bulkVector(i:NX:NX*(NX-2)+i) + bulkVector(NX+i:NX:NX*(NX-1)+i) == 1 );
end

numNoPoints = bulkVector(1) + bulkVector(NX) + bulkVector(end-NX+1) + bulkVector(end) ...
    + sum(bulkVector(1:NX-1) + bulkVector(2:NX) == 2) + sum(bulkVector(end-NX+1:end-1) + bulkVector(end-NX+2:end) == 2) ...
    + sum(bulkVector(1:NX:end-NX) + bulkVector(NX+1:NX:end) == 2) + sum(bulkVector(NX:NX:end-NX) + bulkVector(2*NX:NX:end) == 2);

for i = 0 : NX-2
    numNoPoints = numNoPoints ...
        + sum( bulkVector(i*NX+1:(i+1)*NX-1) + bulkVector(i*NX+2:(i+1)*NX) ...
               + bulkVector((i+1)*NX+1:(i+2)*NX-1) + bulkVector((i+1)*NX+2:(i+2)*NX) == 4 );
end

N_triangles = 2 * N_air;
% N_edges     = 2 * (NX+1)^2 + N_air - numNoEdges;
N_points    = (NX+1)^2 - numNoPoints;

triangles = zeros(N_triangles, 5);
edges = zeros(N_edges, 4);
points = zeros(N_points, 2);

helperPoints = zeros(NX+1,NX+1);

nT = 1;
nE = 1;
nP = 1;

triangles(:,1)  = 3;
triangles(:,2)  = 1;
i = 1;
j = 1;

for  n = 1 : length(bulkVector)
    
    if bulkVector(n)
      j = j + 1;
      if j == NX + 1
        j = 1;
        i = i+1;
      end
      continue
    end
    
    if helperPoints(i,j) == 0
        helperPoints(i,j)   = nP;
        triangles(nT,   3)  = nP;
        triangles(nT+1, 3)  = nP;
        points(nP,:)        = [(j-1)*h, (i-1)*h];
        nP = nP + 1;
    else
        triangles(nT,   3)  = helperPoints(i,j);
        triangles(nT+1, 3)  = helperPoints(i,j);
    end
    
    if helperPoints(i+1,j) == 0
        helperPoints(i+1,j) = nP;
        triangles(nT,   4)  = nP;
        points(nP,:)        = [(j-1)*h, i*h];
        nP = nP + 1;
    else
        triangles(nT,   4)  = helperPoints(i+1,j);
    end
    
    if helperPoints(i,j+1) == 0
        helperPoints(i,j+1) = nP;
        triangles(nT+1, 4)  = nP;
        points(nP,:)        = [j*h, (i-1)*h];
        nP = nP + 1;
    else
        triangles(nT+1, 4)  = helperPoints(i,j+1);
    end
    
    if helperPoints(i+1,j+1) == 0
        helperPoints(i+1,j+1) = nP;
        triangles(nT,   5)  = nP;
        triangles(nT+1, 5)  = nP;
        points(nP,:)        = [j*h, i*h];
        nP = nP + 1;
    else
        triangles(nT,   5)  = helperPoints(i+1,j+1);
        triangles(nT+1, 5)  = helperPoints(i+1,j+1);
    end
    
    if i == 1
        if bulkVector((NX-1)*NX+j)
            edges(nE,:) = [2 2      helperPoints(i,j) helperPoints(i,j+1)];
        else
            edges(nE,:) = [2 888    helperPoints(i,j) helperPoints(i,j+1)];
        end
        nE = nE + 1;
    end
    
    if i == NX
        if bulkVector(j);
            edges(nE,:) = [2 2      helperPoints(i+1,j) helperPoints(i+1,j+1)];
        else
            edges(nE,:) = [2 888    helperPoints(i+1,j) helperPoints(i+1,j+1)];
        end
        nE = nE + 1;
    end
    
    if j == 1
        if bulkVector(i*NX)
            edges(nE,:) = [2 2      helperPoints(i,j) helperPoints(i+1,j)];
        else
            edges(nE,:) = [2 888    helperPoints(i,j) helperPoints(i+1,j)];
        end
        nE = nE + 1;
    end
    
    if j == NX
        if bulkVector((i-1)*NX+1);
            edges(nE,:) = [2 2      helperPoints(i,j+1) helperPoints(i+1,j+1)];
        else
            edges(nE,:) = [2 888    helperPoints(i,j+1) helperPoints(i+1,j+1)];
        end
        nE = nE + 1;
    end
    
    if i ~= 1 && bulkVector(n-NX)
        edges(nE,:)     = [2 2      helperPoints(i,j) helperPoints(i,j+1)];
        nE = nE + 1;
    end
    if i ~= NX && bulkVector(n+NX)
        edges(nE,:)     = [2 2      helperPoints(i+1,j) helperPoints(i+1,j+1)];
        nE = nE + 1;
    end
    if j ~= 1 && bulkVector(n-1)
        edges(nE,:)     = [2 2      helperPoints(i,j) helperPoints(i+1,j)];
        nE = nE + 1;
    end
    if j ~= NX && bulkVector(n+1)
        edges(nE,:)     = [2 2      helperPoints(i,j+1) helperPoints(i+1,j+1)];
        nE = nE + 1;
    end
    
    nT = nT + 2;
    j = j + 1;
    if j == NX + 1
        j = 1;
        i = i+1;
    end
    
end

triangles(:,3:5) = triangles(:,3:5) - 1;
edges(:,3:4) = edges(:,3:4) - 1;
fprintf(file, 'POINTS:');
fprintf(file, '\n%f %f', points.');
fprintf(file, '\nCELLS:');
fprintf(file, '\n%d %d %d %d %d', triangles.');
fprintf(file, '\nFACES:');
fprintf(file, '\n%d %d %d %d', edges.');

end  % function