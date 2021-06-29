clc
close all
clear all


freqAmostragem = 10000;
nPeriodos = 10;
tempoSimbolo_T = 1;
fc = 5 / tempoSimbolo_T;
numBits_k = 2;
numSinais_M = 2^( numBits_k );
indicesAmplitudes_i = randi( numSinais_M, nPeriodos );
amplitude_Pulso = 1;
energiaPulso_Eg = tempoSimbolo_T * amplitude_Pulso;
amplitudeBase_d = sqrt( energiaPulso_Eg / 2 ); 

eixoTemporal = 0:( 1/freqAmostragem ): ( nPeriodos * tempoSimbolo_T ) - ( 1/freqAmostragem );
pulso_g = amplitude_Pulso * ones( freqAmostragem, 1 );
sinalASK_BB = zeros( nPeriodos * freqAmostragem, 1 );

for iPeriodo = 1 : nPeriodos
    
    amplitude = ( 2 * indicesAmplitudes_i( iPeriodo ) - 1 - numSinais_M ) * amplitudeBase_d;
    indiceNoPeriodo = ( iPeriodo - 1 ) * freqAmostragem + 1 : ( iPeriodo ) * freqAmostragem; 
    sinalASK_BB( indiceNoPeriodo ) = amplitude * sqrt( 2 / energiaPulso_Eg ) * amplitude_Pulso;
    
end
figure
subplot(2,1,1);
plot( eixoTemporal, sinalASK_BB, 'Linewidth', 2, 'Color', 'k' );
ylim( [ ( 2 - 1 - numSinais_M - 2 ) * amplitudeBase_d, ( numSinais_M + 1 ) * amplitudeBase_d ] );
grid on;

codigoGray = bin2gray( ( 0 : numSinais_M - 1 )', 'pam', numSinais_M );

propriedadesEixo = gca;
eixoAmplitude = sort( ( 2 * ( 1 : numSinais_M ) - 1 - numSinais_M ) );
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.YTick = eixoAmplitude;
propriedadesEixo.YTickLabel = { strcat( '( ', num2str( dec2base( codigoGray, 2 ) ), ')$\;$',...
    num2str( eixoAmplitude' ) ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 13;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( [ '\textbf{Sinal ASK Banda Base} -- $M =\;$' num2str( numSinais_M ) '$, d = \sqrt{E_g / 2}$' ],...
    'Interpreter', 'Latex', 'FontSize', 15 );
 
sinalASK_PF =  sinalASK_BB .* cos( 2 * pi * fc * eixoTemporal' );

subplot(2,1,2);
plot( eixoTemporal, sinalASK_PF, 'Linewidth', 2, 'Color', 'k' );
ylim( [ ( 2 - 1 - numSinais_M - 2 ) * amplitudeBase_d, ( numSinais_M + 1 ) * amplitudeBase_d ] );
grid on;

propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.YTick = eixoAmplitude;
propriedadesEixo.YTickLabel = { strcat( '( ', num2str( dec2base( codigoGray, 2 ) ), ')$\;$',...
    num2str( eixoAmplitude' ) ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 13;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( [ '\textbf{Sinal ASK Passa Faixa} -- $M =\;$' num2str( numSinais_M ) ', $d = \sqrt{E_g / 2}$' ],...
    'Interpreter', 'Latex', 'FontSize', 15 );

