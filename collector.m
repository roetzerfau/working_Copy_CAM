load sammlung07.mat
display('loaded')

max_runs = 100;
% n = 2;
% m = 30;
%% 
for I = 1:max_runs
fprintf('-------------Schritt %d-------------\n',I);
% display('berechne testdaten \n')
[ews,dist]=dtens_pdist_MainAggregatFolded(0.7);
% display('------------------')
len = size(ews,2);
    for k = 1:len
        [dist,ews,sammlung_dist,sammlung_ews]=sorter(dist,ews,sammlung_dist,sammlung_ews);
    end
    save sammlung07.mat
    display('-------------saved-------------')
end
%%
