function [strucMatrix] = readTxt(filename)

    fileID = fopen(filename,'r');
    formatSpec = '%f';
    strucMatrix = fscanf(fileID,formatSpec);
   
    sizeStruc = sqrt(size(strucMatrix,1));
    strucMatrix = reshape(strucMatrix,[sizeStruc,sizeStruc]);

end