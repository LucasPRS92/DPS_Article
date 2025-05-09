clc
clear all

wk = linspace( -20 , 20 , 1000 );
epsilon = 10;

a = [ 0.9091 0.0178  0.0039 ];
b = [ -1     -0.1087 -0.039 ];
t = [ 1 5 10 ];

[ gw1 ] = evaluateZALMSFunction( wk , epsilon );
[ gw2 ] = evaluatePiecewiseLinearFunction( wk , a , b , t );

set( figure , 'Color' , 'w' )
plot( wk , gw1 , 'r' , 'LineWidth', 2 )
hold on
grid on
plot( wk , gw2 , 'b' , 'LineWidth', 2  )
axis tight


