function [ D ] = assembleGlobInt( g , hatInt , K11DG , K12DG , K21DG , K22DG , ...
    theta1_1 , theta1_2 , theta2_1 , theta2_2 )

[K, N] = size(K11DG);
globInt = cell(4, 1);

for n = 1 : 4
    globInt{n} = sparse(K*N, K*N);
end  % for n

for i = 1 : N
    globInt{1} = globInt{1} + kron( spdiags( K11DG(:,i) , 0 , K , K ) , g.deltaX^2 * hatInt(:,:,i) );
    globInt{2} = globInt{2} + kron( spdiags( K12DG(:,i) , 0 , K , K ) , g.deltaX^2 * hatInt(:,:,i) );
    globInt{3} = globInt{3} + kron( spdiags( K21DG(:,i) , 0 , K , K ) , g.deltaX^2 * hatInt(:,:,i) );
    globInt{4} = globInt{4} + kron( spdiags( K22DG(:,i) , 0 , K , K ) , g.deltaX^2 * hatInt(:,:,i) );
end  % for i

D = zeros(2,2);
D(1,1) = dot( globInt{1} * theta1_1 + globInt{2} * theta1_2 , theta1_1 ) ...
            + dot( globInt{3} * theta1_1 + globInt{4} * theta1_2 , theta1_2 );
D(1,2) = dot( globInt{1} * theta1_1 + globInt{2} * theta1_2 , theta2_1 ) ...
            + dot( globInt{3} * theta1_1 + globInt{4} * theta1_2 , theta2_2 );
D(2,1) = D(1,2);
D(2,2) = dot( globInt{1} * theta2_1 + globInt{2} * theta2_2 , theta2_1 ) ...
            + dot( globInt{3} * theta2_1 + globInt{4} * theta2_2 , theta2_2 );

end  % function