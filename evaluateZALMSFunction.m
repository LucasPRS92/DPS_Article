function [ gw ] = evaluateZALMSFunction( wk , epsilon )

gw = - sign( wk ) ./ ( 1 + epsilon * abs( wk ) );
