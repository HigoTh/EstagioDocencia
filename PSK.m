clc
clear all
close all

freqAmostragem = 10000;
nPeriodos = 10;
tempoSimbolo_T = 1;
fc = 5.5 / tempoSimbolo_T;
numBits_k = 2;
numSinais_M = 2^( numBits_k );
indicesAmplitudes_i = randi( numSinais_M, nPeriodos, 1 );
amplitude_Pulso = 1;
energiaPulso_Eg = tempoSimbolo_T * amplitude_Pulso;
energiaSinal_Es = energiaPulso_Eg / 2;

eixoTemporal = 0:( 1/freqAmostragem ): ( nPeriodos * tempoSimbolo_T ) - ( 1/freqAmostragem );
pulso_g = amplitude_Pulso * ones( freqAmostragem, 1 );
sinalPSK_BBI =  zeros( nPeriodos * freqAmostragem, 1 );
sinalPSK_BBQ =  zeros( nPeriodos * freqAmostragem, 1 );
fase = zeros( nPeriodos, 1 );
faseVec = zeros( nPeriodos * freqAmostragem, 1 );

for iPeriodo = 1 : nPeriodos
    
    fase( iPeriodo ) = ( ( 2 * indicesAmplitudes_i( iPeriodo ) - 1 ) * pi ) / numSinais_M;
    indiceNoPeriodo = ( iPeriodo - 1 ) * freqAmostragem + 1 : ( iPeriodo ) * freqAmostragem; 
    sinalPSK_BBI( indiceNoPeriodo ) = sqrt( ( 2 * energiaSinal_Es ) / energiaPulso_Eg ) *...
        pulso_g * cos( fase( iPeriodo ) );
    sinalPSK_BBQ( indiceNoPeriodo ) = sqrt( ( 2 * energiaSinal_Es ) / energiaPulso_Eg ) *...
        pulso_g * sin( fase( iPeriodo ) );    
    faseVec( indiceNoPeriodo ) = ones( freqAmostragem, 1 ) * fase( iPeriodo );
end

subplot( 3, 2, 1 );
plot( eixoTemporal, sinalPSK_BBI, 'Color', 'k', 'Linewidth', 2 );
ylim( [-2 * amplitude_Pulso, 2 * amplitude_Pulso ] );
grid on

propriedadesEixo = gca;
eixoAmplitude = sort( sqrt( ( 2 * energiaSinal_Es ) / energiaPulso_Eg ) * energiaPulso_Eg * ...
    cos( ( ( 2 * ( 1 : numSinais_M ) - 1 ) * pi ) / numSinais_M ) );
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

for iPeriodo = 1 : nPeriodos
    
    text(  tempoSimbolo_T * ( iPeriodo ) - 0.5 , 1.5 * amplitude_Pulso,...
        [ num2str( rad2deg( fase( iPeriodo ) ) ) '$^{\circ}$' ], 'Interpreter', 'Latex',...
        'HorizontalAlignment','center', 'FontSize', 10 );
    
end

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 13 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 13 );
title( '$s_{I}(t)$', 'Interpreter', 'Latex', 'FontSize', 13 );

subplot( 3, 2, 2 );
plot( eixoTemporal, sinalPSK_BBQ, 'Color', 'k', 'Linewidth', 2 );
ylim( [-2 * amplitude_Pulso, 2 * amplitude_Pulso ] );
grid on

propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

for iPeriodo = 1 : nPeriodos
    
    text(  tempoSimbolo_T * ( iPeriodo ) - 0.5 , 1.5 * amplitude_Pulso,...
        [ num2str( rad2deg( fase( iPeriodo ) ) ) '$^{\circ}$' ], 'Interpreter', 'Latex',...
        'HorizontalAlignment','center', 'FontSize', 10 );
    
end

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 13 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 13 );
title( '$s_{Q}(t)$', 'Interpreter', 'Latex', 'FontSize', 13 );

subplot(3,2,3);
plot( eixoTemporal, sqrt( sinalPSK_BBI.^2 + sinalPSK_BBQ.^2 ) , 'Color', 'k', 'Linewidth', 2 );
ylim( [-2 * amplitude_Pulso, 2 * amplitude_Pulso ] );
grid on;

propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 13 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 13 );
title( '\textbf{Envolt{\''o}ria PSK}', 'Interpreter', 'Latex', 'FontSize', 13 );

subplot(3,2,4);
plot( eixoTemporal, rad2deg( faseVec ) , 'Color', 'k', 'Linewidth', 2 );
ylim( [ min( rad2deg( faseVec ) ) - 10, max( rad2deg( faseVec ) ) + 10 ] );
grid on;

propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.YTick = rad2deg( sort( unique( fase ) ) );
propriedadesEixo.YTickLabel = { strcat( num2str( rad2deg( sort( unique( fase ) ) ) ) ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 13 );
ylabel( 'Fase ($^{\circ}$)', 'Interpreter', 'Latex', 'FontSize', 13 );
title( '\textbf{Fase PSK ($^{\circ}$})', 'Interpreter', 'Latex', 'FontSize', 13 );

subplot(3,2,[5,6]);
sinalPSK_PF = real( sqrt( sinalPSK_BBI.^2 + sinalPSK_BBQ.^2 )' .* ...
    exp( 1i * ( 2 * pi * fc * eixoTemporal + atan( sinalPSK_BBQ ./ sinalPSK_BBI )' ) ) );

plot( eixoTemporal, sinalPSK_PF, 'Color', 'k', 'Linewidth', 2 );
ylim( [ min( rad2deg( faseVec ) ) - 10, max( rad2deg( faseVec ) ) + 10 ] );
grid on;
ylim( [-2.5 * amplitude_Pulso, 2.5 * amplitude_Pulso ] );

propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 13;

codigoGray = bin2gray( ( 0 : numSinais_M - 1 )', 'psk', numSinais_M );
codigoGrayStr = num2str( dec2base( codigoGray, 2 ) );
for iPeriodo = 1 : nPeriodos

    text(  tempoSimbolo_T * ( iPeriodo ) - 0.5 , 1.6 * amplitude_Pulso,...
        [ num2str( rad2deg( fase( iPeriodo ) ) ) '$^{\circ}$' char(10) ...
          '(' codigoGrayStr( indicesAmplitudes_i( iPeriodo ),: ) ')'], 'Interpreter', 'Latex',...
        'HorizontalAlignment','center', 'FontSize', 13 );
    
end

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 13 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 13 );
title( '\textbf{Sinal PSK Passa Faixa} -- $M = 4$, $E_s = E_g / 2$', 'Interpreter', 'Latex', 'FontSize', 13 );