%Calculates the FFT
%Plots the Spectrum
%Plots the PSD
clc;
clear;
set(0, 'DefaultAxesFontSize', 16);

i=1;
j=1;
n=1;
files=dir('/Users/edson/Documents/MATLAB/Sensor_MATLAB/Analise Dados/*');
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
pastas=files(dirFlags);
cd('/Users/edson/Documents/MATLAB/Sensor_MATLAB/Analise Dados/');
delete 'Dataset_Full.CSV'
for i=1:8
        nome_da_pasta = pastas(i).name;
        cd(nome_da_pasta);
        d = dir('*.CSV');
        cd ..;
        for j=1:5
            cd(nome_da_pasta);     
            x_2048=csvread (d(j).name,0,1);
            cd ..;
            k=1;
            for k=1:256:2048
                x=x_2048(k:k+255);
                [L,R] = size(x);
                Fs=4096; %Defined in the mbed code
                %f = Fs*(0:(L/2))/L;
                %f=f(2:(L/2)+1)';
                pw_DC = rms(x)^2;
                x = x - mean(x);
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

    %             pd=makedist('Normal','mu',mean(x),'sigma',desvio);
    %             pdf_x=pdf(pd,x);
    %             pdf_entropy=entropy(pdf_x);
    %             pdf_variance=var(pdf_x);
    %             pdf_kurtosis=kurtosis(pdf_x);
    %             pdf_skewness=skewness(pdf_x);

                y=fft(x,L);
                %y=y(2:(L/2)+1);
                %y_abs=abs(y);
                %y_dB=20*log10(y_abs);
                %plot(f,y_dB);
                xdft = y;
                xdft = xdft(2:L/2+1);
                psdx = (1/(Fs*L)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                psd_mediana=median(psdx);
                psd_kurtosis=kurtosis(psdx);
                psd_skewness=skewness(psdx);
                psd_desvio=std(psdx);
                P2 = 20*log10(abs(y/L));
                P1 = P2(2:L/2+1);
                P1 = 2*P1(1:end-1);
                f = Fs*(1:(L/2-1))/L;
                [M,I]= max(P1);
                beat_frequency=f(I);
                matriz=[pw_AC zcdOut beat_frequency desvio psd_desvio];
                dlmwrite('Dataset_Full.CSV',matriz,'-append');
            end
        end 
end
         
        

%plot no tempo
% tensao = (x.*5.5)./(2^16);
% tempo=(1/Fs:1/Fs:L/Fs);
% 
% subplot(3,1,1);
% plot(tempo,tensao) 
% title('X(t)')
% xlabel('Time (s)')
% ylabel('Voltage (V)')
% 
% subplot(3,1,2);
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('Magnitude (dB)')


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


%velocity=beat_frequency/19.49
%flow_rate=velocity*(0.0025*(pi^2))