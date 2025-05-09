function [ OptimumKappa ] = deriveOptimumKappa2( beta , sigmax2 , w0 , a , b , t , sigmanu2 )

    Delta0 = 1 - beta * sigmax2;
    indices0 =  find( w0 == 0 );
    indicesS =  find( abs( w0 ) < t( end ) );
    indicesS( w0( indicesS ) == 0) = [];
    Q0 = numel( indices0 );
    QS = numel( indicesS );
    b1 = b( 1 );
    [ gw ] = evaluatePiecewiseLinearFunction( w0 , a , b , t );
    gjwSMALL  = gw( indicesS );
    gjwSMALL2 = gjwSMALL .^ 2;

    chi1  = - 2 * Delta0 * b1 * sqrt( 2 / pi );
    chi2  = 2 * beta * sigmax2 * Delta0;
    chi3 = chi2 / ( beta ^ 2 * sigmax2 );
    brevechij = beta * sigmax2 * chi2 ./ ( 2 * Delta0 .* gjwSMALL2 );
    chi4 = Q0 * chi1 / ( 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2 * ( 1 - Q0 * beta ^ 2 * sigmax2 ^ 2 / ( 2 * beta * sigmax2 * Delta0 ) ) );
    Psi1 = Q0 * chi1 * chi3 / ( 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2);
    Omega0 = chi4 * Psi1 / ( Psi1 ); 
    Omega1 = ( - chi1 - 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2 / ( Q0 * chi1 ) * sum( 1 ./ ( brevechij ) ) );
    Omega2 = ( - 4 * beta ^ 3 * sigmax2 ^ 2 * Delta0 * sigmanu2 / chi1 - 8 * beta ^ 4 * sigmax2 ^ 3 * Delta0 ^ 2 * sigmanu2 * QS / ( Q0 * chi1 * chi2 ) );
    Omega3 = 8 * beta ^ 3 * sigmax2 ^ 2 * Delta0 * sigmanu2;
    zeta4 = 0; 
    Omega4 = - Omega0 * Omega1;
    Omega5 = ( chi1 ^ 2 / Omega0 ^ 2 );
    Omega6 = Omega3 / Omega0 ^ 2;

    OptimumKappa = (2^(1/2)*(-(Omega6*(Omega4*(Omega4 + (Omega4^2 - Omega0^4*Omega5)^(1/2)) + Omega0^4*Omega5 - 2*Omega4^2))/(Omega5*(- Omega4^2 + Omega0^4*Omega5)))^(1/2))/2;
 