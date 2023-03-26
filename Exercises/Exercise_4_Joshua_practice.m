clc;close all;clear all

[Y,fs] = audioread('tms.wav');
figure;
pwelch(Y,[],[],[],fs);

SNR = 15;
Y_noise = awgn(Y,SNR);
figure;
pwelch(Y_noise,[],[],[],fs);

fs = 1*fs;
sound(Y,fs)
pause(1)
sound(Y_noise,fs)