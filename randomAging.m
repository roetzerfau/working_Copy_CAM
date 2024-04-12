 function [edgeChargeVector, POMagentAge] = randomAging(POMagentAge, edgeChargeVector)
    
    randHelper = randi(500, size(edgeChargeVector));
    randHelper = randHelper > 499;
    
    edgeChargeVector((POMagentAge > 50) & randHelper) = 0;
 end