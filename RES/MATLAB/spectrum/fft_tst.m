clear;
close all;

dt = 0.01;
T = 1000;
N = round(T/dt);
t = 0:dt:(T-dt);

y = [ones(1, 1000) zeros(1, N-1000)];
y = sin(2*t) + sin(t*sqrt(2)) + 3*cos(15*t) + 2;
[A, w] = fft_m(y, dt);
plot(w, A);


% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 3500;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% plot_time = 150;
% 
% S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
% X = S + 0.5*randn(size(t));
% 
% getFig('t (milliseconds)', 'X(t)', 'Signal Corrupted with Zero-Mean Random Noise');
% plot(1000*t(1:plot_time),X(1:plot_time));
% 
% Y = fft(X);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% 
% getFig('f (Hz)','|P1(f)|','Single-Sided Amplitude Spectrum of X(t)');
% plot(f,P1);
% %plot(Y/L);
% 
% Y = fft(S);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% getFig('f (Hz)','|P1(f)|','Single-Sided Amplitude Spectrum of S(t)');
% plot(f,P1);
