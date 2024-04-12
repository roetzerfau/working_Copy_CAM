function printInfoTUM(k,bulkVector,POMconcVector, concPOMAgent, edgeChargeVector, POMsolidEdgeList,...
    numFreePOMparticles, numEdgeTypes, totalPOMinputConc, totalPOMoutputConc, sumExcessPOM, POMagentInput,...
    POMocclusion_total, POMocclusion_attractive)

if k == 0
   flag = 'w';
else
   flag = 'a';
end

fileName    = 'txtdata/POMconc.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(POMconcVector));
fclose(fileID_k);

fileName    = 'txtdata/POMconcNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(POMconcVector > 0));
fclose(fileID_k);

fileName    = 'txtdata/concPOMAgent.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(concPOMAgent));
fclose(fileID_k);

fileName    = 'txtdata/concPOMAgentNNZ.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(concPOMAgent > 0));
fclose(fileID_k);

fileName    = 'txtdata/attractiveEdges.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sum(edgeChargeVector));
fclose(fileID_k);

fileName    = 'txtdata/edgeTypes.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e %e %e %e %e %e \n', k, numEdgeTypes(1), numEdgeTypes(2), ...
    numEdgeTypes(3), numEdgeTypes(4), numEdgeTypes(5), numEdgeTypes(6));
fclose(fileID_k);

surfaceCoverage = 0;
for i = 1 : length(POMsolidEdgeList)
    surfaceCoverage = surfaceCoverage + length(POMsolidEdgeList{i});
end
fileName    = 'txtdata/surfaceCoverage.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, surfaceCoverage);
fclose(fileID_k);

fileName    = 'txtdata/POMocclusion_total.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, POMocclusion_total);
fclose(fileID_k);

fileName    = 'txtdata/POMocclusion_attractive.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, POMocclusion_attractive);
fclose(fileID_k);

fileName    = 'txtdata/numFreePOMparticles.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, numFreePOMparticles);
fclose(fileID_k);

fileName    = 'txtdata/totalPOMinput.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, totalPOMinputConc);
fclose(fileID_k);

fileName    = 'txtdata/totalPOMoutput.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, totalPOMoutputConc);
fclose(fileID_k);

fileName    = 'txtdata/sumExcessPOM.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, sumExcessPOM);
fclose(fileID_k);

fileName    = 'txtdata/POMagentInput.txt';
fileID_k = fopen(fileName,flag);
fprintf(fileID_k, '%d %e \n', k, POMagentInput);
fclose(fileID_k);

end


