
clc;
clear;
set(0, 'DefaultAxesFontSize', 16);
%x=csvread ('/Users/edson/Documents/MATLAB/Sensor de adubo/SAMPLES_NoFlow.CSV',0,1);
x=csvread('/Users/edson/Documents/MATLAB/Sensor de adubo/Analise Dados/85g_s/SAMPLES 5.CSV',0,1);
x=x-mean(x);
pw_DC = rms(x)^2;
[L,R] = size(x);
Fs=4096; %Defined in the mbed code
pw_AC = rms(x)^2;
%s = spectrogram(x);
%spectrogram(x,[],[],[],Fs,'yaxis');
zcd = dsp.ZeroCrossingDetector;
zcdOut = step(zcd,x);
x_normalized = x/max(abs(x));
kt=kurtosis(x);
e = entropy(x_normalized);
variance=var(x);
mediana=median(x);
desvio=std(x);
y=fft(x,L);
P2 = 20*log10(abs(y/L));
P1 = P2(2:L/2+1);
P1 = 2*P1(1:end-1);
f = Fs*(1:(L/2-1))/L;
[M,I]= max(P1);
beat_frequency=f(I);


%plot no tempo
tensao = (x.*5.5)./(2^16);
tempo=(1/Fs:1/Fs:L/Fs);

% subplot(3,1,1);
plot(tempo,tensao) 
axis([0 0.5 -6.0 4.0])
title('X(t)')
xlabel('Time (s)')
ylabel('Voltage (V)')

% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('Magnitude (dB)')
% 
% 
% xdft = y;
% xdft = xdft(2:L/2+1);
% psdx = (1/(Fs*L)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 1:Fs/L:Fs/2;
% subplot(3,1,3);
% plot(freq,10*log10(psdx))
% grid on
% title('PSD Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')
