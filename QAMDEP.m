clc
clear all
close all

T = 1;
n_period = 4;
sig_a = 1;

f = linspace( 0, 4 * ( 1 / T ), 10000 );

S = T * sig_a^2 * ( sin( pi * f * T ) ./ ( pi * f * T ) ).^2;

plot( f, 10*log10( S ), 'Color', 'k', 'Linewidth', 2 );

ylim( [ -50, max( 10*log10( S ) ) ] );
grid on;

xlabel( 'Frequ{\^{e}}ncia (Hz)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( '$S_{\tilde{v}}(f)$ (dB)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( '\textbf{DEP do Sinal $\tilde{v}(t)$}', 'Interpreter', 'Latex', 'FontSize', 15 );

propriedadesEixo = gca;
propriedadesEixo.XTick = 1 : n_period;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), '/T' ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 13;