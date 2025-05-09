function [ aj , bj ] = returnAjBj( wk , a , b , t )

P = numel( a );
N = numel( wk );

aj = zeros( size( wk ) );
bj = zeros( size( wk ) );

indices = cell( P , 1 );

for p = 1 : P

    if ( p == 1 )

        indices{ 1 } = find( ( abs( wk ) >= 0 ) & ( abs( wk ) < t( 1 ) ) );
        aj( indices{ 1 } ) = a( 1 );
        bj( indices{ 1 } ) = b( 1 ) * sign( wk( indices{ 1 } ) ); 

    else

        indices{ p } = find( ( abs( wk ) >= t( p - 1 ) ) & ( abs( wk ) < t( p ) ) );
        aj( indices{ p } ) = a( p );
        bj( indices{ p } ) = b( p ) * sign( wk( indices{ p } ) ); 

    end

end


