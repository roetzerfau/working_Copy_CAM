function [particleList, particleSizesHelper] = createPOMParticleListNew (POMVector) 

NX = sqrt(length(POMVector));
POMVectorFlag  = -POMVector;
particleCounter = 0;
waitList = [];

POMInd = find(POMVector == 1)';

for bulk = POMInd % iterate over all bulk cells
    if POMVectorFlag( bulk ) < 0
        particleCounter = particleCounter + 1;
        waitList = [waitList, bulk];        
        while isempty(waitList) == 0
            currentBulk = waitList(1);
            waitList = waitList(2:end);
            if POMVectorFlag(currentBulk) < 0
                POMVectorFlag(currentBulk) = particleCounter;
            end            
            possibleNeighbors = stencil( NX, NX, currentBulk, 1); 
            possibleNeighbors = possibleNeighbors(2:end);
            neighbors = possibleNeighbors(POMVector(possibleNeighbors) == 1);
            for neighborBulk = neighbors
                if POMVectorFlag( neighborBulk ) < 0
                    POMVectorFlag( neighborBulk ) = particleCounter;
                    waitList = [waitList, neighborBulk];
                end
            end
        end
    end
end


particleList = cell(1,particleCounter);
for particle =  1 : particleCounter
    particleList{ particle } = sort(find(POMVectorFlag == particle)'); 
end
particleSizesHelper = 0;
end
