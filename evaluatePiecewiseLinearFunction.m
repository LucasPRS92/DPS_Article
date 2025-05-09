function [ gw ] = evaluatePiecewiseLinearFunction( wk , a , b , t )

P = numel( a );
N = numel( wk );

gw = zeros( size ( wk ) );

indices = cell( P , 1 );

for p = 1 : P

    if ( p == 1 )

        indices{ 1 } = find( ( abs( wk ) >= 0 ) & ( abs( wk ) < t( 1 ) ) );
        gw( indices{ 1 } ) = a( 1 ) * wk( indices{ 1 } ) + b( 1 ) * sign( wk( indices{ 1 } ) ); 

    else

        indices{ p } = find( ( abs( wk ) >= t( p - 1 ) ) & ( abs( wk ) < t( p ) ) );
        gw( indices{ p } ) = a( p ) * wk( indices{ p } ) + b( p ) * sign( wk( indices{ p } ) ); 

    end

end





