 function POMagentAge = calculatePOMagentAge(parameters, POMagentAge, edgeChargeVector, concPOMAgent )
    POMagentAge(edgeChargeVector > 0) = POMagentAge(edgeChargeVector > 0) + 1;
    POMagentAge(concPOMAgent > parameters.POMagentMin) = 0;
 end