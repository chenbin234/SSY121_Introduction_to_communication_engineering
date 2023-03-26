% An example showing how to load and play a waveform,
% Plot Power spectral density using Welch's periodogram averaging method
clc;  close all; clear all

[Y, fs] = audioread('tms.wav'); %read in the signal
% Plot and look at the signal Y
figure; pwelch(Y,[] ,[] ,[],fs) % help pwelch, Welch's periodogram averaging

% Add noise
SNR = 25;  % signal to noise ratio in dB
Y_noise = awgn(Y, SNR);          % try different SNRs
figure; 
pwelch(Y_noise,[] ,[] ,[],fs) 

%% Play the waveforms
fs = 1*fs;
sound(Y,fs)            % Is fs necessary? Try different fs
pause(1)
sound(Y_noise,fs)

% Generate your own waveforms 