function [ MSD ] = runZALMS( w0 , numberOfIterations , sigmanu2 , beta , kappa , epsilon , numberOfRepeats );

N = numel( w0 );
MSD = zeros( numberOfIterations , 1 );

for repeat = 1 : numberOfRepeats

    wk = zeros( N , 1 );

    x = randn( N + numberOfIterations - 1 , 1 );
    d = filter( w0 , 1 , x );
    d = d + sqrt( sigmanu2 ) * randn( size( d ) );

    for k = N : numberOfIterations + N - 1

        xk = x( k: -1 : k - N + 1 );
        yk = transpose( wk ) * xk;
        ek = d( k ) - yk;

        wk = wk + beta * xk * ek + kappa * evaluateZALMSFunction( wk , epsilon );

        MSD( k - N + 1 ) = MSD( k - N + 1 ) + norm( wk - w0 ) ^ 2 / numberOfRepeats;

    end

end

