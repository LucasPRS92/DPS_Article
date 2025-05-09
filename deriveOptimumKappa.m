function [ OptimumKappa ] = deriveOptimumKappa( beta , sigmax2 , w0 , a , b , t , sigmanu2 )
    kappa = 0;

    Delta0 = 1 - beta * sigmax2;
    indices0 =  find( w0 == 0 );
    indicesS =  find( abs( w0 ) < t( end ) );
    indicesS( w0( indicesS ) == 0) = [];
    indicesL = find( abs( w0 ) > t( end ) );
    Q0 = numel( indices0 );
    QS = numel( indicesS );
    QL = numel( indicesL );
    lado1 = 0;
    lado2 = 0;
    a1 = a( 1 );
    b1 = b( 1 );
    [ gw ] = evaluatePiecewiseLinearFunction( w0 , a , b , t );
    gjwSMALL  = gw( indicesS );
    gjwSMALL2 = gjwSMALL .^ 2;
    [ aj , ~ ] = returnAjBj( w0( indicesS ) , a , b , t );


            eta10 = 2 * beta * sigmax2 * Delta0;


            chi1  = - 2 * Delta0 * b1 * sqrt( 2 / pi );
            eta20 = chi1 * kappa;
            eta30 = -beta^2*sigmanu2*sigmax2                     ;
        
            chi2  = 2 * beta * sigmax2 * Delta0;
            overlinechij = - aj .^ 2;
            eta3i = chi2                                                  + overlinechij * kappa ^ 2;

            chi3 = chi2 / ( beta ^ 2 * sigmax2 );
            hatchij = overlinechij / ( beta ^ 2 * sigmax2 );
            eta1i = 1 ./ ( chi3 + hatchij * kappa ^ 2);
            eta1i = beta ^ 2 * sigmax2 ^ 2 ./ eta3i;

            brevechij = beta * sigmax2 * chi2 ./ ( 2 * Delta0 .* gjwSMALL2 );
            dotchij   = beta * sigmax2 * overlinechij ./ ( 2 * Delta0 .* gjwSMALL2 ); 
            eta2i = kappa ^ 2 ./ ( brevechij  + dotchij * kappa ^ 2 );
            chi4 = Q0 * chi1 / ( 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2 * ( 1 - Q0 * beta ^ 2 * sigmax2 ^ 2 / ( 2 * beta * sigmax2 * Delta0 ) ) );
            Psi1 = Q0 * chi1 * chi3 / ( 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2);
            dotPsij = Q0 * chi1 * hatchij / ( 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2 );
            Omega0 = chi4 * Psi1 / ( Psi1 - QS * chi4 );
            zeta1 = 1 / ( Omega0 * kappa );
                                                  
 zeta2 = ( - chi1 - 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2 / ( Q0 * chi1 ) * sum( 1 ./ ( brevechij ) ) ) * kappa + ( - 4 * beta ^ 3 * sigmax2 ^ 2 * Delta0 * sigmanu2 / chi1 - 8 * beta ^ 4 * sigmax2 ^ 3 * Delta0 ^ 2 * sigmanu2 * QS / ( Q0 * chi1 * chi2 ) ) * ( 1 / kappa );
            Omega1 = ( - chi1 - 8 * beta ^ 2 * sigmax2 ^ 2 * Delta0 ^ 2 / ( Q0 * chi1 ) * sum( 1 ./ ( brevechij ) ) );
            Omega2 = ( - 4 * beta ^ 3 * sigmax2 ^ 2 * Delta0 * sigmanu2 / chi1 - 8 * beta ^ 4 * sigmax2 ^ 3 * Delta0 ^ 2 * sigmanu2 * QS / ( Q0 * chi1 * chi2 ) );
            zeta2 = Omega1 * kappa + Omega2 * ( 1 / kappa );
  


            Omega3 = 8 * beta ^ 3 * sigmax2 ^ 2 * Delta0 * sigmanu2;
            zeta3 = chi1 ^ 2 * kappa ^ 2 + Omega3;
            zeta4 = 4 * beta ^ 2 * eta10 * sigmax2 ^ 2;
            
            Omega4 = ( zeta4 * Omega0 ^ 2 / 2 - Omega0 * Omega1 );
            Omega5 = ( chi1^2/Omega0^2 - Omega1 / Omega0 * zeta4 );
            Omega6 = ( Omega3 / Omega0 ^ 2 - Omega2 * zeta4 / Omega0 );

OptimumKappa = (2^(1/2)*(-(Omega6*(Omega4*(Omega4 + (Omega4^2 - Omega0^4*Omega5)^(1/2)) + Omega0^4*Omega5 - 2*Omega4^2))/(Omega5*(- Omega4^2 + Omega0^4*Omega5)))^(1/2))/2;
 