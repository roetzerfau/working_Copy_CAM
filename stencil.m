function [ stencil ] = stencil( NX , NZd , candidates , layers )
%STENCIL Code to create a stencil of an arbitrary number of layers

%% Error handling

if any( candidates( : ) < 0 ) || any( candidates( : ) == 0 )
    error('Stencil is only applicalble to real positive integers')
end

if layers <= 0 || mod( layers , 1 ) ~= 0
    error('layer has to be a real positive integer')
end

%% Stencil building

% compute the number of elements of a stencil with "layers" layers
stencil          = zeros( length( candidates ) , 2 * ( layers ) * ( layers + 1 ) + 1 );
% the first entry of a stencil is allways the cell where it is centered
% stencil( : , 1 ) = candidates;
% stencil = stencil_layer_1(NX,NZd,candidates);
stencil(:,1:5) =  stencil_layer_1(NX,NZd,candidates);   
if layers > 1
    for m = 1 : layers-1
        sten = stencil_layer_1(NX,NZd,stencil(:,2*(m-1)*m+2));
        stencil(:,2*m*(m+1)+2) = sten(:,2);
        stencil(:,2*m*(m+1)+3) = sten(:,3);
        stencil(:,2*m*(m+1)+4) = sten(:,4);
        
        for i = 1 : 2*(m+1)-3
            sten1 = stencil_layer_1(NX,NZd,stencil(:,2*(m-1)*m+2+2*i-1));
            sten2 = stencil_layer_1(NX,NZd,stencil(:,2*(m-1)*m+2+2*i));
            stencil(:,2*m*(m+1)+4+2*i-1)=sten1(:,3);
            stencil(:,2*m*(m+1)+4+2*i)=sten2(:,4);
        end
        sten = stencil_layer_1(NX,NZd,stencil(:,2*m*(m+1)+1));
        stencil(:,2*(m+1)*(m+2)-1) = sten(:,3);
        stencil(:,2*(m+1)*(m+2)  ) = sten(:,4);
        stencil(:,2*(m+1)*(m+2)+1) = sten(:,5);

    end
end
end


function [stencil] = stencil_layer_1( NX , NZd , candidates)
layers = 1;
stencil          = zeros( length( candidates ) , 2 * ( layers ) * ( layers + 1 ) + 1 );
% the first entry of a stencil is allways the cell where it is centered
stencil( : , 1 ) = candidates;
numberOfElements = 1;

for i = 1 : layers
    for j = -i : i

        if abs( j ) == i
            stencil( : , numberOfElements + 1 ) = NX * mod( floor( ( candidates + ( j * NX ) - 1 ) / NX ) , NZd ) + mod( candidates - NX - 1 , NX ) + 1;
            numberOfElements = numberOfElements + 1;
        else
            leftIdx  = - ( i - abs( j ) ) - 1;
            rightIdx = + ( i - abs( j ) );

            stencil( : , numberOfElements + 1 ) = NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates + leftIdx, NX ) + 1 ) + ( j * NX ) - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates + leftIdx , NX ) + 1 ) + ( j * NX ) - 1 , NX ) + 1;
            stencil( : , numberOfElements + 2 ) = NX * mod( floor( ( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + rightIdx ) + ( j * NX ) - 1 ) / NX ) , NZd ) + mod( ( NX * floor( ( ( candidates - 1 ) / NX ) ) + mod( candidates , NX ) + rightIdx ) + ( j * NX ) - 1 , NX ) + 1;
            
            numberOfElements = numberOfElements + 2;
        end % if 
    end % for j
end % for i 


end