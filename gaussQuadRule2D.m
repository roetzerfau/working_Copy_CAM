function [ Q1 , Q2 , W ] = gaussQuadRule2D( maxExp )

% Note: maxExp is the maximum exponent of x and y!

switch maxExp
    case {0 , 1}
        Q1 = 0.5 ;
        Q2 = 0.5 ;
        W  = 1 ;
    case {2 , 3}
        Q1 = [ 0.211324865405187 ; 0.211324865405187 ; 0.788675134594813 ; ... 
               0.788675134594813 ];
        Q2 = [ 0.211324865405187 ; 0.788675134594813 ; 0.211324865405187 ; ...
               0.788675134594813 ];
        W  = 0.25 * ones(1 , 4);
    case {4 , 5}
        Q1 = [ 0.112701665379258 ; 0.112701665379258 ; 0.112701665379258 ; ... 
               0.5               ; 0.5               ; 0.5               ; ...
               0.887298334620742 ; 0.887298334620742 ; 0.887298334620742 ]; 
        Q2 = [ 0.112701665379258 ; 0.5               ; 0.887298334620742 ; ...
               0.112701665379258 ; 0.5               ; 0.887298334620742 ; ...
               0.112701665379258 ; 0.5               ; 0.887298334620742 ];
        W  = [ 0.077160493827160 , 0.123456790123457 , 0.077160493827160 , ...
               0.123456790123457 , 0.197530864197531 , 0.123456790123457 , ...
               0.077160493827160 , 0.123456790123457 , 0.077160493827160 ];
    case {6 , 7}
        Q1 = [ 0.069431844202974 ; 0.069431844202974 ; 0.069431844202974 ; ...
               0.069431844202974 ; 0.330009478207572 ; 0.330009478207572 ; ...
               0.330009478207572 ; 0.330009478207572 ; 0.669990521792428 ; ...
               0.669990521792428 ; 0.669990521792428 ; 0.669990521792428 ; ...
               0.930568155797026 ; 0.930568155797026 ; 0.930568155797026 ; ...
               0.930568155797026 ]; 
        Q2 = [ 0.069431844202974 ; 0.330009478207572 ; 0.669990521792428 ; ...
               0.930568155797026 ; 0.069431844202974 ; 0.330009478207572 ; ...
               0.669990521792428 ; 0.930568155797026 ; 0.069431844202974 ; ...
               0.330009478207572 ; 0.669990521792428 ; 0.930568155797026 ; ...
               0.069431844202974 ; 0.330009478207572 ; 0.669990521792428 ; ...
               0.930568155797026 ];
        W  = [ 0.030250748321401 , 0.056712962962963 , 0.056712962962963 , ...
               0.030250748321401 , 0.056712962962963 , 0.106323325752674 , ...
               0.106323325752674 , 0.056712962962963 , 0.056712962962963 , ...
               0.106323325752674 , 0.106323325752674 , 0.056712962962963 , ...
               0.030250748321401 , 0.056712962962963 , 0.056712962962963 , ...
               0.030250748321401 ];
    case {8 , 9}
        Q1 = [ 0.046910077030668 ; 0.046910077030668 ; 0.046910077030668 ; ...
               0.046910077030668 ; 0.046910077030668 ; 0.230765344947159 ; ...
               0.230765344947159 ; 0.230765344947159 ; 0.230765344947159 ; ...
               0.230765344947159 ; 0.5               ; 0.5               ; ...
               0.5               ; 0.5               ; 0.5               ; ...
               0.769234655052842 ; 0.769234655052842 ; 0.769234655052842 ; ...
               0.769234655052842 ; 0.769234655052842 ; 0.953089922969332 ; ...
               0.953089922969332 ; 0.953089922969332 ; 0.953089922969332 ; ...
               0.953089922969332 ]; 
        Q2 = [ 0.046910077030668 ; 0.230765344947159 ; 0.5               ; ...
               0.769234655052842 ; 0.953089922969332 ; 0.046910077030668 ; ...
               0.230765344947159 ; 0.5               ; 0.769234655052842 ; ...
               0.953089922969332 ; 0.046910077030668 ; 0.230765344947159 ; ...
               0.5               ; 0.769234655052842 ; 0.953089922969332 ; ...
               0.046910077030668 ; 0.230765344947159 ; 0.5               ; ...
               0.769234655052842 ; 0.953089922969332 ; 0.046910077030668 ; ...
               0.230765344947159 ; 0.5               ; 0.769234655052842 ; ...
               0.953089922969332 ];
        W  = [ 0.014033587215607 , 0.028350000000000 , 0.033696268096880 , ...
               0.028350000000000 , 0.014033587215607 , 0.028350000000000 , ...
               0.057271351055998 , 0.068071633137688 , 0.057271351055998 , ...
               0.028350000000000 , 0.033696268096880 , 0.068071633137688 , ...
               0.080908641975309 , 0.068071633137688 , 0.033696268096880 , ...
               0.028350000000000 , 0.057271351055998 , 0.068071633137688 , ...
               0.057271351055998 , 0.028350000000000 , 0.014033587215607 , ...
               0.028350000000000 , 0.033696268096880 , 0.028350000000000 , ...
               0.014033587215607 ];
end  % switch

end  % function