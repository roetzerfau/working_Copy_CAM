% Does not work for correctly, not correctly implemented for particle_index
% > 4. Only a problem if values on edges, e.g. charge vary from horizontal
% to horizontal edge, or vertical to vertical edge.
function bulk_trafo = get_particle_trafo(candidate, g, NZd, rotationDirection, particle_index)
bulkSize = length(candidate);
if particle_index == 3
    sten = stencil(g.NX,NZd,candidate(1),1);
    if ~ismember(sten(4),candidate) % vertikal
        if rotationDirection == 1 
            bulk_trafo = [4 1 8 7 2 11 3 10 5 6 12 9];
        else
            bulk_trafo = [3 7 6 1 10 6 4 2 12 8 5 11];
        end
    else                            % horizontal
        if rotationDirection == 1 
            bulk_trafo = [4 8 1 7 11 3 2 10 6 5 12 9];
        else
            bulk_trafo = [2 5 7 1 9 10 4 3 12 8 6 11];
        end
    end 
elseif particle_index == 1
    sten = stencil(g.NX,NZd,candidate(1),1);
    if ~ismember(sten(4),candidate) % vertikal
        if rotationDirection == 1 
            bulk_trafo = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
        else
            bulk_trafo = [1 3 2 5 4 7 6 9 8 11 10 13 12 15 14 17 16];
        end
    else                            % horizontal
        if rotationDirection == 1 
            bulk_trafo = [1 3 2 5 4 7 6 9 8 11 10 13 12 15 14 17 16];
        else
            bulk_trafo = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
        end
    end 
elseif particle_index == 4
    sten = stencil(g.NX,NZd,candidate(1),1);
    if ~ismember(sten(4),candidate) % vertikal
        if rotationDirection == 1 
            bulk_trafo = [4 1 8 7 2 12 3 11 5 16 6 15 9 20 10 19 13 24 14 23 17 28 18 27 21 32 22 31 25 36 26 35 29 40 30 39 31 44 34 43 37 48 38 47 41 52 42 51 45 56 46 55 49 59 50 58 53 54 60 57];
        else
            bulk_trafo = [3 7 6 1 11 10 4 2 15 14 8 5 19 18 12 9 23 22 16 13 27 26 20 17 31 30 24 21 35 34 28 25 39 38 32 29 43 42 36 33 47 46 40 37 51 50 44 41 55 54 48 45 58 57 52 49 60 56 53 59];
        end
    else                            % horizontal
        if rotationDirection == 1 
            bulk_trafo = [2 5 7 1 9 11 4 3 13 15 8 6 17 19 12 10 21 23 16 14 25 27 20 18 29 31 24 22 33 35 28 26 37 39 32 30 41 43 36 34 45 47 40 38 49 51 44 42 53 55 48 46 57 58 52 50 60 56 54 59];
        else
            bulk_trafo = [4 8 1 7 12 3 2 11 16 6 5 15 20 10 9 19 24 14 13 23 28 18 17 21 32 22 21 31 36 26 25 35 40 30 29 39 44 34 33 43 48 38 37 47 52 42 41 51 56 46 45 55 59 50 49 58 54 53 60 57];
        end
    end 
elseif particle_index == 2
    sten = stencil(g.NX,NZd,candidate(1),1);
    if ~ismember(sten(4),candidate) % vertikal
        if rotationDirection == 1 
            bulk_trafo = [4 1 8 7 2 12 3 11 5 16 6 15 9 20 10 19 13 24 14 23 17 28 18 27 21 32 22 31 25 36 26 35 29 40 30 39 31 44 34 43 37 48 38 47 41 52 42 51 45 56 46 55 49 60 50 59 53 64 54 63 57 67 58 66 61 62 68 65];
        else
            bulk_trafo = [3 7 6 1 11 10 4 2 15 14 8 5 19 18 12 9 23 22 16 13 27 26 20 17 31 30 24 21 35 34 28 25 39 38 32 29 43 42 36 33 47 46 40 37 51 50 44 41 55 54 48 45 59 58 52 49 63 62 56 53 66 65 60 57 68 64 61 67];
        end
    else                            % horizontal
        if rotationDirection == 1 
            bulk_trafo = [2 5 7 1 9 11 4 3 13 15 8 6 17 19 12 10 21 23 16 14 25 27 20 18 29 31 24 22 33 35 28 26 37 39 32 30 41 43 36 34 45 47 40 38 49 51 44 42 53 55 48 46 57 59 52 50 61 63 56 54 65 66 60 58 68 64 62 67];
        else
            bulk_trafo = [4 8 1 7 12 3 2 11 16 6 5 15 20 10 9 19 24 14 13 23 28 18 17 21 32 22 21 31 36 26 25 35 40 30 29 39 44 34 33 43 48 38 37 47 52 42 41 51 56 46 45 55 60 50 49 59 64 54 53 63 67 58 57 66 62 61 68 65];
        end
    end 
elseif particle_index == 6
    bulk_trafo = [1];
elseif particle_index == 7
    bulk_trafo = [1 2 3 4 5 6 7 8 9 10 11 12];   
elseif particle_index == 8 
    bulk_trafo = 1:36; 
elseif particle_index == 9 || particle_index == 10 || particle_index == 11 || particle_index == 12  
    bulk_trafo = 1:30;   
elseif particle_index == 13 || particle_index == 14 || particle_index == 15 || particle_index == 16  
    bulk_trafo = 1:60;    
elseif particle_index == 17 || particle_index == 18 || particle_index == 19 || particle_index == 20  
    bulk_trafo = 1:48;        
elseif bulkSize == 800
    sten = stencil(g.NX,NZd,candidate(1),4);
    if ~ismember(sten(34),candidate) % vertikal
        bulk_trafo = [3 5 6 1 2 4];
    else                            % horizontal
        bulk_trafo = [1 3 4 2 6 5];
    end  
else
    error('wrong candidate');
end
end