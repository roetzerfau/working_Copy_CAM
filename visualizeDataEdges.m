function visualizeDataEdges( g , repLagr , varName, fileName, tLvl, edgeType )

[K, N] = size(repLagr);
%% open file
fileName    = ['vtk/' , fileName, '.', num2str(tLvl), '.vtu'];
file        = fopen(fileName, 'wt');
%% header
fprintf(file, '<?xml version="1.0"?>\n');
fprintf(file, '<VTKFile type="UnstructuredGrid" versig.on="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(file, '  <UnstructuredGrid>\n');
%% points & cells
if edgeType == 1
    P1 = g.coordV(reshape(g.V0E', 2*length(g.V0E),1),1);
    P2 = g.coordV(reshape(g.V0E', 2*length(g.V0E),1),2);
elseif edgeType == 2
    P1 = g.coordV(reshape(g.V0CE', 2*length(g.V0CE),1),1);
    P2 = g.coordV(reshape(g.V0CE', 2*length(g.V0CE),1),2);
end
numP = 2;
id = 3;

fprintf(file, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',K*numP,K);
fprintf(file, '      <Points>\n');
fprintf(file, '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(file, '          %.3e %.3e %.3e\n',  [P1, P2, zeros(numP*K, 1)]');
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Points>\n');
fprintf(file, '      <Cells>\n');
fprintf(file, '        <DataArray type="Int32" Name="connectivity" format="ascii">\n');
fprintf(file,'           '); fprintf(file,'%d ', 0:K*numP-1);
fprintf(file, '\n        </DataArray>\n');
fprintf(file, '        <DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(file,'           %d\n', numP:numP:numP*K);
fprintf(file, '        </DataArray>\n');
fprintf(file, '        <DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(file,'           %d\n', id*ones(K, 1));
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </Cells>\n');
%% data

dataLagr = kron(repLagr, [1;1])';

fprintf(file, '      <PointData Scalars="%s">\n', varName);
fprintf(file, '        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n', varName);
fprintf(file, '          %.3e\n', dataLagr);
fprintf(file, '        </DataArray>\n');
fprintf(file, '      </PointData>\n');
%% footer
fprintf(file, '    </Piece>\n');
fprintf(file, '  </UnstructuredGrid>\n');
fprintf(file, '</VTKFile>\n');
%% close file
fclose(file);
disp(['Data written to ' fileName])
end

