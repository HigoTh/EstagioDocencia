clc
clear all
close all

Fs = 10000;
nPeriodos = 10;
T = 1;
fc = 5.5 / T;
k = 4;
M = 2^( k );
iI = randi( sqrt( M ), nPeriodos, 1 );
iQ = randi( sqrt( M ), nPeriodos, 1 );
ampPulso = 1;
Eg = T * ampPulso;
d = sqrt( Eg / 2 );
index = ( iI - 1 ) * sqrt( M ) + ( sqrt( M ) - iQ ) + 1;

tempo = 0:( 1/Fs ): ( nPeriodos * T ) - ( 1/Fs );

pulso_g = ampPulso * ones( Fs, 1 );
sinalQAM_BBI =  zeros( nPeriodos * Fs, 1 );
sinalQAM_BBQ =  zeros( nPeriodos * Fs, 1 );
fasePeriodo = zeros( nPeriodos, 1 );
envoltPeriodo = zeros( nPeriodos, 1 );

for iPeriodo = 1 : nPeriodos
    
    A = ( 2 * iI( iPeriodo ) - sqrt( M ) - 1 ) * d;
    B = ( 2 * iQ( iPeriodo ) - sqrt( M ) - 1 ) * d;
    indPeriodo = ( iPeriodo - 1 ) * Fs + 1 : ( iPeriodo ) * Fs;
    
    sinalQAM_BBI( indPeriodo ) = A * sqrt( 2 / Eg ) * pulso_g;
    sinalQAM_BBQ( indPeriodo ) = B * sqrt( 2 / Eg ) * pulso_g;

    fasePeriodo( iPeriodo ) = atan2( B, A );
    envoltPeriodo( iPeriodo ) = sqrt( A^2 + B^2 ) * sqrt( 2 / Eg );
        
end

figure
subplot(3,2,1);
plot( tempo, sinalQAM_BBI, 'Linewidth', 2, 'Color', 'k' );
ylim( [ ( 2 - 1 - sqrt(M) - 2 ) * d, ( sqrt(M) + 1 ) * d ] );
grid on;

codigoGray = bin2gray( ( 0 : sqrt( M ) - 1 )', 'qam', sqrt( M ) );

propriedadesEixo = gca;
eixoAmplitude = ( 2 * ( 1 : sqrt( M ) ) - sqrt( M ) - 1 );
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.YTick = eixoAmplitude;
propriedadesEixo.YTickLabel = { strcat( '( ', num2str( dec2base( codigoGray, 2 ) ), ')$\;$',...
    num2str( eixoAmplitude' ) ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( [ '\textbf{Sinal QAM banda base (em fase)} -- $M =\;$' num2str( M ) ],...
    'Interpreter', 'Latex', 'FontSize', 12 );

subplot(3,2,2);
plot( tempo, sinalQAM_BBQ, 'Linewidth', 2, 'Color', 'k' );
ylim( [ ( 2 - 1 - sqrt(M) - 2 ) * d, ( sqrt(M) + 1 ) * d ] );
grid on;

propriedadesEixo = gca;
eixoAmplitude = ( 2 * ( 1 : sqrt( M ) ) - sqrt( M ) - 1 );
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.YTick = eixoAmplitude;
propriedadesEixo.YTickLabel = { strcat( '( ', num2str( dec2base( codigoGray, 2 ) ), ')$\;$',...
    num2str( eixoAmplitude' ) ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( [ '\textbf{Sinal QAM banda base (em quadratura)} -- $M =\;$' num2str( M )  ],...
    'Interpreter', 'Latex', 'FontSize', 12 );

eixoFase = zeros( M, 1 );
eixoEnvolt = zeros( M, 1 );
k = 0;
for i = 1 : sqrt( M )
    for j = 1 : sqrt( M )
        k = k + 1;
        Ai = ( 2 * i - sqrt( M ) - 1 ) * d;
        Bi = ( 2 * j - sqrt( M ) - 1 ) * d;
        eixoFase( k ) = atan2( Bi , Ai );
        eixoEnvolt( k ) = sqrt( Ai^2 + Bi^2 ) * sqrt( 2 / Eg );
    end
end

envoltoria = sqrt( sinalQAM_BBI.^2 + sinalQAM_BBQ.^2 );
subplot(3,2,3);
plot( tempo, envoltoria, 'Linewidth', 2, 'Color', 'k' );
ylim( [ min( sqrt( 2 * eixoAmplitude.^2 ) ) - sqrt(2) * d, max( sqrt( eixoAmplitude.^2 ) ) + 2 * sqrt(2) * d ] );
grid on;
propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
propriedadesEixo.YTick = sort( unique( eixoEnvolt ) );
propriedadesEixo.YTickLabel = { num2str( propriedadesEixo.YTick', 3 ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( [ '\textbf{Envolt{\''o}ria QAM} -- $M =\;$' num2str( M ) ],...
    'Interpreter', 'Latex', 'FontSize', 12 );


fase = atan2( sinalQAM_BBQ , sinalQAM_BBI );
subplot(3,2,4);
plot( tempo, rad2deg( fase ), 'Linewidth', 2, 'Color', 'k' );
ylim( [ min( rad2deg( eixoFase ) ) - 10, max( rad2deg( eixoFase ) ) + 10 ] );
grid on;
propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.YTick = rad2deg( sort( unique( eixoFase ) ) );
propriedadesEixo.YTickLabel = { strcat( num2str( rad2deg( sort( unique( eixoFase ) ) ), 4 ) ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 9;

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Fase ($^{\circ}$)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( [ '\textbf{Fase QAM} -- $M =\;$' num2str( M ) ],...
    'Interpreter', 'Latex', 'FontSize', 12 );

sinalQAM = envoltoria .* cos( 2* pi * fc * tempo' + fase );
subplot(3,2,[5,6]);
plot( tempo, sinalQAM, 'Linewidth', 2, 'Color', 'k' );
hold on
plot( tempo, envoltoria, '--', 'Linewidth', 2, 'Color', 'r' );
hold on
plot( tempo, -envoltoria, '--', 'Linewidth', 2, 'Color', 'r' );
ylim( [ -( max( envoltoria ) + 3 * d ), 2.5 * max( envoltoria ) ] );

propriedadesEixo = gca;
propriedadesEixo.XTick = 0 : nPeriodos;
propriedadesEixo.XTickLabel = { strcat( num2str( propriedadesEixo.XTick' ), 'T' ) };
propriedadesEixo.XTickLabel( 1 ) = { num2str( 0 ) };
% propriedadesEixo.YTick = sort( unique( eixoEnvolt ) );
% propriedadesEixo.YTickLabel = { num2str( propriedadesEixo.YTick', 2 ) };
propriedadesEixo.TickLabelInterpreter = 'Latex';
propriedadesEixo.FontSize = 12;

codigoGray = bin2gray( ( 0 : M - 1 )', 'qam', M );
grayStr = num2str( dec2base( codigoGray, 2 ) );

for iPeriodo = 1 : nPeriodos
    
    txt = [ '$V_i =\;$' num2str( envoltPeriodo( iPeriodo ), 3 ) ' V' char(10) '$\phi_i =\;$'...
        num2str( rad2deg( fasePeriodo( iPeriodo ) ), 4 ) '$\;^{\circ}$' ...
        char(10) '(' grayStr( index( iPeriodo ), : ) ')' ];
    
    text(  T * ( iPeriodo ) - 0.5 , 2.4 * max( envoltoria ),txt,...
        'Interpreter', 'Latex','HorizontalAlignment','center', 'VerticalAlignment', 'top', 'FontSize', 9 );
    
end

xlabel( 'Tempo (s)', 'Interpreter', 'Latex', 'FontSize', 15 );
ylabel( 'Amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 15 );
title( [ '\textbf{Envolt{\''o}ria QAM} -- $M =\;$' num2str( M ) ],...
    'Interpreter', 'Latex', 'FontSize', 12 );

close all
x = (0:M - 1)';
y = bin2gray(x,'qam',M);
hMod = comm.RectangularQAMModulator;
symbols = constellation(hMod);
h = scatterplot(symbols,1,0,'s');
title( '\textbf{Varia{\c{c}}{\~a}o Temporal dos S{\''i}mbolos na Constela{\c{c}}{\~a}o QAM}',...
    'Interpreter','Latex');
xlabel('');
ylabel('');
for k = 1:16
    text(real(symbols(k))-0.3,imag(symbols(k))+0.3,...
        dec2base(y(k),2,4));
end
axis([-4 4 -4 4]);
hold on
% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for iPeriodo = 1 : nPeriodos-1
    % Draw plot for y = x.^n
    plot(real(symbols(index(iPeriodo:iPeriodo + 1))),imag(symbols(index(iPeriodo:iPeriodo + 1))),'Color','r','LineWidth',2); 
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if iPeriodo == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end



