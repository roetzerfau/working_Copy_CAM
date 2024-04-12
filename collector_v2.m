
%%
load sammlung.mat
display('sammlung loaded')
aufruf = aufruf+1; %% Anzahl, wie oft collector_v2 schon aufgerufen wurde
max_runs = 1000;
pfad = ['Verteilungen/Startverteilungen/Aufruf',num2str(aufruf)];
%% 
mkdir(pfad)
for I = 1:max_runs
fprintf('-------------Schritt %d von %d-(0.8)--------\n',I,max_runs);
[ews,dist]=dtens_pdist_MainAggregatFolded(0.8);
copyfile('solu.0.vtu',[pfad,'/aufruf_',num2str(aufruf),'_run_',num2str(I),'.vtu']);
display(['Startverteilung solu.0.vtu moved to ',pfad,'/aufruf_',num2str(aufruf),'_run_',num2str(I),'.vtu'])

len = size(ews,2);
    for k = 1:len
        [dist,ews,sammlung_dist,sammlung_ews]=sorter_v2(dist,ews,sammlung_dist,sammlung_ews,I,k,aufruf);
    end
    save sammlung.mat
    display('-------------saved-------------')
end
%%
