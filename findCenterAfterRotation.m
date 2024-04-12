function centerInd = findCenterAfterRotation(center, N, height, width, rotationDirection)
% Input: 
% center = index of center cell of particle
% height = height of particle
% width = width of particle
% N = domain of size N * N
% rotationDirection = index of direction in which the particle rotates
% Output:
% centerInd = index of center cell after rotation of particle
    
% indices start from 1 at left lower corner to N*N in right upper corner, 
% going row by row

% range of distLeftBoundary, distRightBoundary : 0 to N - 1
    distLeftBoundary = mod(center - 1,N);
    distRightBoundary = mod(N*N - center,N);
    
    
    if width > height % horizontal particle
        
        % number of cells the center cell moves left/right AND up/down,
        % e.g. for a horizontal particle of size 2x6, the center cell moves
        % exactly 2 cells in two directions
        rotDist = (width - height)/2;
        
        % four directions for particle to rotate
        switch rotationDirection 
        
        % rotation left + up
        case 1
            % center cell moves over left boundary
            if rotDist > distLeftBoundary
                centerInd = center + N*floor(rotDist) + (N-ceil(rotDist)); 
            % center cell does not cross boundary    
            else
                centerInd = center + N*floor(rotDist) - ceil(rotDist);
            end
        % rotation right + up    
        case 2
            % center cell moves over right boundary
            if rotDist > distRightBoundary
                centerInd = center + N*floor(rotDist) - (N-floor(rotDist));
            % center cell does not cross boundary     
            else
                centerInd = center + N*floor(rotDist) + floor(rotDist);
            end
        % rotation left + down     
        case 3
        % center cell moves over left boundary    
            if rotDist > distLeftBoundary
                centerInd = center - N*ceil(rotDist) + (N-ceil(rotDist)); 
            % center cell does not cross boundary     
            else
                centerInd = center - N*ceil(rotDist) - ceil(rotDist);
            end
        % rotation right + down 
        case 4
            % center cell moves over right boundary
            if rotDist > distRightBoundary
                centerInd = center - N*ceil(rotDist) - (N-floor(rotDist)); 
            % center cell does not cross boundary     
            else
                centerInd = center - N*ceil(rotDist) + floor(rotDist);
            end
        end
        
    else % vertical particle
        
        % analogeously rotation distance for vertically placed particle
        rotDist = (height - width)/2;
        
        switch rotationDirection 
        case 1
            if rotDist > distLeftBoundary
                centerInd = center + N*ceil(rotDist) + (N-floor(rotDist)); 
            else
                centerInd = center + N*ceil(rotDist) - floor(rotDist);
            end
        case 2
            if rotDist > distRightBoundary
                centerInd = center + N*ceil(rotDist) - (N-ceil(rotDist)); 
            else
                centerInd = center + N*ceil(rotDist) + ceil(rotDist);
            end
        case 3           
            if rotDist > distLeftBoundary
                centerInd = center - N*floor(rotDist) + (N-floor(rotDist)); 
            else
                centerInd = center - N*floor(rotDist) - floor(rotDist);
            end

        case 4
            if rotDist > distRightBoundary
                centerInd = center - N*floor(rotDist) - (N-ceil(rotDist)); 
            else
                centerInd = center - N*floor(rotDist) + ceil(rotDist);
            end
        end
    end
    
    
    % center cell moved below lower boundary, due to periodicity is moved
    % to top of domain
    if centerInd <= 0
        centerInd = centerInd + N*N;
    end
    
    % center cell moved above upper boundary
    if centerInd > N*N
        centerInd = centerInd - N*N;
    end
    
end
