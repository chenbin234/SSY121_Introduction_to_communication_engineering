% a)
close all;
clc;
clear all;

flag =1; %0 = PN sequence, 1 = Barker, 2 = given sequence
SNR = 20;
switch flag
    case 0        
        preamble = [1 1 1 -1 1 -1 -1];   % PN sequence
        
    case 1
        preamble = [1 1 1 -1 -1 1 -1];     % Barker code, https://en.wikipedia.org/wiki/Barker_code
        
    case 2
        preamble = [1 -1 1 -1 1 -1 1];   % given sequence
    case 3
        preamble = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; % 13 length Barker code
end

plot(xcorr(preamble))

Tp = length(preamble);          % length of preamble
delay = randi([1,20], 1);       % Generate a random delay between 1 and 20
rx = awgn([zeros(1, delay) preamble], SNR, 'measured');     % insert the sequence after the delay
corr = conv(rx, fliplr(conj(preamble)));   % correlate the sequence and received vector
figure; plot(corr, '.-r')       % plot correlation

[tmp, Tmax] = max(corr);         %find location of max correlation
Tx_hat = Tmax - length(preamble);  %find delay (it will take length(pream
fprintf('delay = %d, estimation = %d\n', delay, Tx_hat)


% b
preamble_error = preamble;
preamble_error(1) = -preamble_error(1);
preamble_error(7) = -preamble_error(7);

rx = awgn([zeros(1, delay) preamble_error], SNR, 'measured');     % insert the sequence after the delay
corr = conv(rx, fliplr(preamble));   % correlate the sequence and received vector
figure; plot(corr, '.-r')       % plot correlation

[tmp, Tmax] = max(corr);         %find location of max correlation
Tx_hat = Tmax - length(preamble_error);  %find delay
fprintf('delay = %d, estimation = %d\n', delay, Tx_hat)



%% 
% Is the choice of the sequence important for estimating the delay?
% Try different sequences and see if it answers the question above. 
% Try, plot(xcorr(seq)) on different sequences and observe the shape and
% width besides the point of maximum correlation coefficient
