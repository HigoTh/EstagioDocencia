clc
clear all
close all

MVec = 256;
EbNoVec = linspace( -40, 60, 10000 ) + 10 * log10( log2( MVec ) );
BW = 200e6;
% MVec = 2.^( 2 : 8 );
X = [ 29.0, 22.8, 25.8, 13.8, 16.8, 4.8, 7.8 ];
Y = [ 7.98, 7.11, 7.76, 4.35, 5.18, 1.69, 1.94 ];

for i = 1 : size( MVec, 2 )

    berTheory = berawgn( EbNoVec,'qam',MVec );
    H = -( 1 - berTheory ) .* log2( 1 - berTheory ) - ( berTheory .* log2( berTheory ) );
    C = log2( MVec ) * ( 1 - H );
    plot( EbNoVec, C, 'LineWidth', 2 );
    hold on
    
end
grid on
ax = gca;
ax.TickLabelInterpreter = 'Latex';
ax.FontSize = 11;
% plot( X,Y, 'o' );
% ylim( [ 1e-6, 1] );
% xlim( [ -10, 50 ] );

% MStr = num2str( MVec' );
% legendString = { strcat( '$$M =\;$$', MStr ) };
% legend( legendString, 'Interpreter', 'Latex' );
% 
% title( '\textbf{Probabilidade de Erro de bit -- QAM}', 'Interpreter', 'Latex', 'FontSize', 15  );
% xlabel( '$E_{b}/N_0\;$(dB)', 'Interpreter', 'Latex', 'FontSize', 13  );
% ylabel( '$P(e)$', 'Interpreter', 'Latex', 'FontSize', 13  );


% QAM = [ 9.9, 14.0, 18.6 ];
% PSK = [ 9.9, 18.1, 28.3 ];
% ASK = [ 13.8, 23.1, 33.5 ];
% 
% M = [ 4, 16, 64 ];
% Rb_W = log2( M ) / 2;
% 
% semilogy( QAM, Rb_W, '--s' );
% hold on
% semilogy( ASK, Rb_W, '--s' );
% hold on
% semilogy( PSK, Rb_W, '--s' );
% hold on



