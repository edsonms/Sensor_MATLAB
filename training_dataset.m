%Calculates the FFT
%Plots the Spectrum
%Plots the PSD
clc;
clear;
set(0, 'DefaultAxesFontSize', 16);

files=dir('/Users/edson/Documents/MATLAB/Sensor_MATLAB/Analise Dados/*');
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
pastas=files(dirFlags);
cd('/Users/edson/Documents/MATLAB/Sensor_MATLAB/Analise Dados/');
delete 'dataset_training.CSV'
for i=1:8
        nome_da_pasta = pastas(i).name;
        cd(nome_da_pasta);
        d = dir('*.CSV');
        cd ..;
        for j=1:3
            cd(nome_da_pasta);     
            x_2048=csvread (d(j).name,0,1);
            cd ..;
            k=1;
            for k=1:256:2048
                x=x_2048(k:k+255);
                [L,R] = size(x);
                Fs=4096; %Defined in the mbed code

                pw_DC = rms(x)^2;
                x = x - mean(x);
                pw_AC = rms(x)^2;
                pw_AC = pw_AC/57017000;

                zcd = dsp.ZeroCrossingDetector;
                zcdOut_int = step(zcd,x);
                zcdOut = cast(zcdOut_int,'double');
                x_normalized = x/max(abs(x));
                kt=kurtosis(x);
                e = entropy(x_normalized);
                variance=var(x);
                mediana=median(x);
                desvio=std(x);

                y=fft(x,L);

                xdft = y;
                xdft = xdft(2:L/2+1);
                psdx = (1/(Fs*L)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);

                psd_desvio=std(psdx);
                P2 = 20*log10(abs(y/L));
                P1 = P2(2:L/2+1);
                P1 = 2*P1(1:end-1);
                f = Fs*(1:(L/2-1))/L;
                [M,I]= max(P1);
                beat_frequency=f(I);
                %matriz=[0.5 pw_AC/57017000 zcdOut/596 beat_frequency/108 desvio/7553 psd_desvio/149030 str2num(nome_da_pasta)/100];
                zcdOut=zcdOut/596;
                beat_frequency=beat_frequency/108;
                desvio = desvio/7553;
                psd_desvio=psd_desvio/149030;
                frequencia_esperada=str2num(nome_da_pasta)/100;
                
                matriz=[pw_AC zcdOut beat_frequency desvio psd_desvio frequencia_esperada];
                dlmwrite('dataset_training.CSV',matriz,'precision',10,'-append');
            end
        end 
end
         
i=1;
j=4;
files=dir('/Users/edson/Documents/MATLAB/Sensor_MATLAB/Analise Dados/*');
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
pastas=files(dirFlags);
cd('/Users/edson/Documents/MATLAB/Sensor_MATLAB/Analise Dados/');
delete 'dataset_testing.CSV'
for i=1:8
        nome_da_pasta = pastas(i).name;
        cd(nome_da_pasta);
        d = dir('*.CSV');
        cd ..;
        for j=4:5
            cd(nome_da_pasta);     
            x_2048=csvread (d(j).name,0,1);
            cd ..;
            k=1;
            for k=1:256:2048
                x=x_2048(k:k+255);
                [L,R] = size(x);
                Fs=4096; %Defined in the mbed code

                pw_DC = rms(x)^2;
                x = x - mean(x);
                pw_AC = rms(x)^2;
                pw_AC = pw_AC/57017000;

                zcd = dsp.ZeroCrossingDetector;
                zcdOut_int = step(zcd,x);
                zcdOut = cast(zcdOut_int,'double');
                x_normalized = x/max(abs(x));
                kt=kurtosis(x);
                e = entropy(x_normalized);
                variance=var(x);
                mediana=median(x);
                desvio=std(x);

                y=fft(x,L);

                xdft = y;
                xdft = xdft(2:L/2+1);
                psdx = (1/(Fs*L)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);

                psd_desvio=std(psdx);
                P2 = 20*log10(abs(y/L));
                P1 = P2(2:L/2+1);
                P1 = 2*P1(1:end-1);
                f = Fs*(1:(L/2-1))/L;
                [M,I]= max(P1);
                beat_frequency=f(I);
                zcdOut=zcdOut/596;
                beat_frequency=beat_frequency/108;
                desvio = desvio/7553;
                psd_desvio=psd_desvio/149030;
                frequencia_esperada=str2num(nome_da_pasta)/100;
                matriz=[pw_AC zcdOut beat_frequency desvio psd_desvio frequencia_esperada];
                dlmwrite('dataset_testing.CSV',matriz,'precision',10,'-append');
            end
        end 
end        
