clear all; close all;
[Tfilename, path]=uigetfile('.txt','Open'); 
if(path==0)
  return
end
path_=strcat(path,Tfilename);
rawData = importdata(path_);

%========================================
%                Load data
%========================================

emg1_data = rawData(:,3);
emg2_data = rawData(:,6);
emg3_data = rawData(:,9);

emg1_data = emg1_data*5/1023;
emg2_data = emg2_data*5/1023;
emg3_data = emg3_data*5/1023;

Fs = 1000;

emg1_data(~isfinite(emg1_data))=0;
emg2_data(~isfinite(emg2_data))=0;
emg3_data(~isfinite(emg3_data))=0;

emg1_data = nonzeros(emg1_data);
emg2_data = nonzeros(emg2_data);
emg3_data = nonzeros(emg3_data);

Ty1 = (0:numel(emg1_data)-1)/Fs;
Ty2 = (0:numel(emg2_data)-1)/Fs;
Ty3 = (0:numel(emg3_data)-1)/Fs;

%========================================
%          low, high pass filter
%========================================

% X1 = emg1_data; % 미분된신호 재정의
% n1 = 5;  % 차수
% Wn1 = 10; % 주파수
% Fn = Fs/2; % 샘플링데이터/2
% ftype = 'low'; % type=low
% [b,a]=butter(n1, Wn1/Fn, ftype);
% emg1_lowpass = filter(b, a, X1);
% 
% X2 = emg1_lowpass; % low-pass신호 재정의
% Wn3 = 420; % 주파수
% ftype1 = 'high'; % type=high
% [e,f] = butter(n1, Wn3/Fn, ftype1);
% emg1_low_high_pass = filter(e, f, X2);
% 
% X3 = emg2_data;
% emg2_lowpass = filter(b, a, X3);
% 
% X4 = emg2_lowpass;
% emg2_low_high_pass = filter(e, f, X4);
% 
% X5 = emg3_data;
% emg3_lowpass = filter(b, a, X5);
% 
% X6 = emg3_lowpass;
% emg3_low_high_pass = filter(e, f, X6);
% 
% figure(1)
% subplot(3,1,1)
% plot(Ty1, emg1_low_high_pass, 'b')
% xlabel('Time[s]')
% ylabel('Data')
% 
% subplot(3,1,2)
% plot(Ty2, emg2_low_high_pass, 'b')
% xlabel('Time[s]')
% ylabel('Data')
% 
% subplot(3,1,3)
% plot(Ty3, emg3_low_high_pass, 'b')
% xlabel('Time[s]')
% ylabel('Data')

%========================================
%                   emg
%========================================

% norm_emg1 = emg1_low_high_pass - min_norm_emg;
% norm_emg1 = norm_emg1 / (max_norm_emg - min_norm_emg);
% 
% norm_emg2 = emg2_low_high_pass - min_norm_emg;
% norm_emg2 = norm_emg2 / (max_norm_emg - min_norm_emg);
% 
% norm_emg3 = emg3_low_high_pass - min_norm_emg;
% norm_emg3 = norm_emg3 / (max_norm_emg - min_norm_emg);

%========================================
%                   FFT
%========================================

% a = emg1_data((Fs*0+1):Fs*10);
% b = emg1_data((Fs*10+1):Fs*20);

emg11_data = emg1_data(1:2000);
emg22_data = emg2_data(1:2000);
emg33_data = emg3_data(1:2000);

emg4_data = emg1_data(1001:3000);
emg5_data = emg2_data(2001:4000);
emg6_data = emg3_data(14001:16000);

L1 = length(emg11_data);
L2 = length(emg22_data);
L3 = length(emg33_data);

emg1_FFT = fft(emg11_data);
emg2_FFT = fft(emg22_data);
emg3_FFT = fft(emg33_data);

emg4_FFT = fft(emg4_data);
emg5_FFT = fft(emg5_data);
emg6_FFT = fft(emg6_data);

emg1_P2 = abs(emg1_FFT/L1);
emg1_P1 = emg1_P2(1:L1/2+1);
emg1_P1(2:end-1) = 2*emg1_P1(2:end-1);

emg2_P4 = abs(emg2_FFT/L2);
emg2_P3 = emg2_P4(1:L2/2+1);
emg2_P3(2:end-1) = 2*emg2_P3(2:end-1);

emg3_P6 = abs(emg3_FFT/L3);
emg3_P5 = emg3_P6(1:L3/2+1);
emg3_P5(2:end-1) = 2*emg3_P5(2:end-1);

f1 = Fs*(0:(L1/2))/L1;
f2 = Fs*(0:(L2/2))/L2;
f3 = Fs*(0:(L3/2))/L3;

ff1 = (0:L1-1)*(Fs/L1); 
ff2 = (0:L2-1)*(Fs/L2); 
ff3 = (0:L3-1)*(Fs/L3); 

emg1_power = abs(emg1_FFT).^2/L1;
emg2_power = abs(emg2_FFT).^2/L2;
emg3_power = abs(emg3_FFT).^2/L3;

emg4_power = abs(emg4_FFT).^2/L1;
emg5_power = abs(emg5_FFT).^2/L2;
emg6_power = abs(emg6_FFT).^2/L3;

% med_freq1 = medfreq(emg1_P1, Fs);
% med_freq2 = medfreq(emg1_power, Fs);
% med_freq3 = medfreq(emg3_FFT);

figure(2)
subplot(3,1,1)
plot(ff1,emg4_power) 
% title('Single-Sided Amplitude Spectrum of X(t)')
xlim([1 500])
ylim([0 1])
xlabel('f (Hz)')
ylabel('power spectrum')

subplot(3,1,2)
plot(ff2,emg5_power) 
xlim([1 500])
ylim([0 1])
% title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('power spectrum')

subplot(3,1,3)
plot(ff3,emg6_power) 
xlim([1 500])
ylim([0 1])
% title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('power spectrum')

cums1 = cumsum(emg1_power);
cums2 = cumsum(emg2_power);
cums3 = cumsum(emg3_power);

med_freq1 = medfreq(cums1, ff1);
med_freq2 = medfreq(cums2, ff2);
med_freq3 = medfreq(cums3, ff3);

% figure(3)
% subplot(3,1,1)
% plot(ff1, cums1);
% subplot(3,1,2)
% plot(ff2, cums2);
% subplot(3,1,3)
% plot(ff3, cums3);

% figure(3)
% % subplot(2,1,1)
% plot(ff3,emg3_power) 
% xlim([1 500])
% ylim([0 1])
% xlabel('F (Hz)')
% ylabel('Power Spectrum')
% 
% subplot(2,1,2)
% plot(ff3,emg6_power) 
% xlim([1 500])
% ylim([0 1])
% xlabel('F (Hz)')
% ylabel('Power Spectrum')