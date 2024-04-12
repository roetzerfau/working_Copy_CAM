function [particleList, particleContent] = particleInfoNew(bulkVector,particle1List,particle2List,particle3List,particle4List,particle5List, ...
    particle6List,fluid)
% function based on BFS
if nargin == 7
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
            possibleNeighbors = stencil( NX, NX, currentBulk, 1); possibleNeighbors = possibleNeighbors(2:end);
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
if fluid == 1
    particleContent = 0;
%     particlePerimeter = 0;
%     particleRadius = 0;
    return
end

particleContent = cell(1,particleCounter);
for particle = 1 : particleCounter
    tmp_particleInd = particleList{ particle };
    while ~isempty(tmp_particleInd)
        
        [row,~] = find(particle1List == tmp_particleInd(1) ); %check particle1List
        flag = 1;
        if isempty(row)
            [row,~] = find(particle2List == tmp_particleInd(1) ); %check particle2List
            flag = 2;
        end
        if isempty(row)
            [row,~] = find(particle3List == tmp_particleInd(1) ); %check particle3List
            flag = 3;
        end
        if isempty(row)
            [row,~] = find(particle4List == tmp_particleInd(1) ); %check particle4List
            flag = 4;
        end
        if isempty(row)
            [row,~] = find(particle5List == tmp_particleInd(1) ); %check particle5List
            flag = 5;
        end
        if isempty(row)
            [row,~] = find(particle6List == tmp_particleInd(1) ); %check particle6List
            flag = 6;
        end

        assert(~isempty(row),'Teilchen nicht gefunden')
        assert(length(row)==1,'Zu viele Teilchen gefunden')
        particleContent{ particle } = [particleContent{particle};[flag,row]];
        switch flag
            case 1 %Quartz
                tmp_particleInd = setdiff(tmp_particleInd,particle1List(row,:),'stable');
            case 2 %Goethit
                tmp_particleInd = setdiff(tmp_particleInd,particle2List(row,:),'stable');
            case 3 %Illite
                tmp_particleInd = setdiff(tmp_particleInd,particle3List(row,:),'stable');
            case 4
                tmp_particleInd = setdiff(tmp_particleInd,particle4List(row,:),'stable');
            case 5
                tmp_particleInd = setdiff(tmp_particleInd,particle5List(row,:),'stable');
            case 6
                tmp_particleInd = setdiff(tmp_particleInd,particle6List(row,:),'stable'); 
        end
    end
end
% 
% %TESTING
% summ = 0;
% for particle = 1: length(particleList)
%     summ = summ + length(particleList{particle});
% end
% if summ ~= sum(bulkVector)
%     fprintf('es wurden nicht alle Teilchen entdeckt');
% end
% for particle = 1: length(particleList)
%     test_particle = [];
%     for subparticle = 1: size(particleContent{particle},1)
%         type = particleContent{particle}(subparticle,1);
%         switch type
%             case 1
%                 test_particle = [test_particle, QuartzList(particleContent{particle}(subparticle,2),:)]; 
%             case 2
%                 test_particle = [test_particle, GoethitList(particleContent{particle}(subparticle,2),:)]; 
%             case 3
%                 test_particle = [test_particle, IlliteList(particleContent{particle}(subparticle,2),:)]; 
%         end
%                
%     end
%     test_particle = sort(test_particle);
%     if ~isequal(test_particle,particleList{particle})
%         fprintf('particleContent stimmt nicht');
%     end
% end
end