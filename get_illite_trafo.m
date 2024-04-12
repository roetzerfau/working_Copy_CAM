function bulk_trafo = get_illite_trafo(candidate, g, NZd)
bulkSize = length(candidate);
if bulkSize == 6
    sten = stencil(g.NX,NZd,candidate(1),1);
    if ~ismember(sten(3),candidate) % vertikal
        bulk_trafo = [3 5 6 1 2 4];
    else                            % horizontal
        bulk_trafo = [1 3 4 2 6 5];
    end 
elseif bulkSize == 24
    % noch implementieren
elseif bulkSize == 96
    % noch implementieren
else
    error('wrong candidate');
end
end