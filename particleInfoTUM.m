function [particleList, particleContent] = particleInfoTUM(bulkVector, solidParticleList, POMParticleList, fluid)
% function based on BFS
if nargin == 3
    fluid = 0;
end
particlePerimeter = zeros(1,length(bulkVector));


NX = sqrt(length(bulkVector));
bulkVectorFlag  = -bulkVector;
particleCounter = 0;
waitList = [];

bulkInd = find(bulkVector == 1)';

for bulk = bulkInd % iterate over all bulk cells
    if bulkVectorFlag( bulk ) < 0
        particleCounter = particleCounter + 1;
        waitList = [waitList, bulk];        
        while isempty(waitList) == 0
            currentBulk = waitList(1);
            waitList = waitList(2:end);
            if bulkVectorFlag(currentBulk) < 0
                bulkVectorFlag(currentBulk) = particleCounter;
            end            
            possibleNeighbors = stencil( NX, NX, currentBulk, 1); 
            possibleNeighbors = possibleNeighbors(2:end);
            neighbors = possibleNeighbors(bulkVector(possibleNeighbors) == 1);
            particlePerimeter( particleCounter ) = particlePerimeter( particleCounter ) + 4 - length(neighbors);
            for neighborBulk = neighbors
                if bulkVectorFlag( neighborBulk ) < 0
                    bulkVectorFlag( neighborBulk ) = particleCounter;
                    waitList = [waitList, neighborBulk];
                end
            end
        end
    end
end



particleList = cell(1,particleCounter);
for particle =  1 : particleCounter
    particleList{ particle } = sort(find(bulkVectorFlag == particle)'); 
end

% particleContent = 0;
% return

% Problematic, since solidParticleList is a cell array, differently than
% particle1List etc.

particleContent = cell(1,particleCounter);
for particle = 1 : particleCounter
    tmp_particleInd = particleList{ particle };
    while ~isempty(tmp_particleInd)
        
        particleFound = 0;
        
        for i = 1 : length(solidParticleList)
           rowHelper = find(solidParticleList{i} ==tmp_particleInd(1) );
           if (~isempty(rowHelper))
               row = i;
               flag = 1;
               particleFound = 1;
           end
        end
        
        if particleFound ~= 1
            for i = 1 : length(POMParticleList)
               rowHelper = find(POMParticleList{i} ==tmp_particleInd(1) );
               if (~isempty(rowHelper))
                   row = i;
                   flag = 2;
                   particleFound = 1;
               end
            end
        end


        assert(~isempty(row),'Teilchen nicht gefunden')
        assert(length(row)==1,'Zu viele Teilchen gefunden')
        particleContent{ particle } = [particleContent{particle};[flag,row]];
        switch flag
            case 1 
                tmp_particleInd = setdiff(tmp_particleInd,solidParticleList{row},'stable');
            case 2 
                tmp_particleInd = setdiff(tmp_particleInd,POMParticleList{row},'stable');
        end
    end
end

end