function PorosityDependency(  )

POR = (1:-0.01:0.5);

realPor     = zeros(size(POR));
lambda_1_0  = zeros(size(POR));
lambda_2_0  = zeros(size(POR));
lambda_1_50 = zeros(size(POR));
lambda_2_50 = zeros(size(POR));

i = 1;

while i <= length(POR)
    [real_porosity, K1, K2] = MainAggregatFoldedHelper( POR(i) );
    if sum(sum(isnan(K1))) > 0
        continue;
    elseif sum(sum(isnan(K2))) > 0
        continue;
    end
    a0  = eigs(K1);
    a50 = eigs(K2);
    if a0(1) > 1.1 || a0(2)>1.1 || a50(1)>1.1 || a50(2) > 1.1
        continue;
    end
    if a0(1) > a50(1) + 0.01 || a0(2) > a50(2) + 0.01
        continue;
    end
    realPor(i)      = real_porosity;
    lambda_1_0(i)   = a0(1);
    lambda_2_0(i)   = a0(2);
    lambda_1_50(i)  = a50(1);
    lambda_2_50(i)  = a50(2);
    i = i + 1;
end

fid = fopen('lambda_1^0', 'wt' );
fprintf( fid, '%18.14f %18.14f\n', [realPor.', lambda_1_0.'].' );
fclose( fid );
fid = fopen('lambda_2^0', 'wt' );
fprintf( fid, '%18.14f %18.14f\n', [realPor.', lambda_2_0.'].' );
fclose( fid );
fid = fopen('lambda_1^50', 'wt' );
fprintf( fid, '%18.14f %18.14f\n', [realPor.', lambda_1_50.'].' );
fclose( fid );
fid = fopen('lambda_2^50', 'wt' );
fprintf( fid, '%18.14f %18.14f\n', [realPor.', lambda_2_50.'].' );
fclose( fid );

end