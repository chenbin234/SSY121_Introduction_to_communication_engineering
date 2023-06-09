% Signal space and constellation
% Ex 5
clc; clear all; close all

% Input papameters
% For passband implemetation sampling frequency is a necessary parameter
fs = 100; % sampling frequency [Hz]
Tsamp = 1/fs;
Rb = 100; % bit rate [bit/sec]
N = 100; % number of bits to transmit
% Constellation or bit to symbol mapping
const = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % Constellation 1 - QPSK/4-QAM
symbs = N/log2(length(const));

pulse_flag = 2; % 0 = rect, 1 = rc, 2 = rrc

% s = exp(1i*((0:3)*pi/2 + pi/4)); % Constellation 1 - same constellation generated as PSK
scatterplot(const); grid on;                            % Constellation visualization

M = length(const);                                      % Number of symbols in the constellation
bpsymb = log2(M);                                        % Number of bits per symbol
fsymb = Rb/bpsymb;                                          % Symbol rate [symb/s]
Tsymb = 1/fsymb;
fsfd = fs/fsymb;                                       % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity) [samples/symb]
a = randsrc(1,N,[0 1]);                             % Information bits
m = buffer(a, bpsymb)';                           % Group bits into bits per symbol
m = bi2de(m, 'left-msb')'+1;           % Bits to symbol index
x = const(m);                                     % Look up symbols using the indices

x_upsample = upsample(x, fsfd);                     % Space the symbols fsfd apart, to enable pulse shaping using conv.

% Now we have a signal ready to be shaped.

switch pulse_flag
    case 0 %Rectangular
        span = 0.5;
        pulse = ones(1, fsfd);                   % Generate a rectangular pulse
        pulse = pulse/norm(pulse);
        s = conv(pulse,x_upsample);    % Each symbol replaced by the pulse shape and added, this is our pulse train.
        
    case 1 %RC
        span = 6;
        [pulse, t] = rcpuls(0.35,1/fsymb,fs,span);
        s = conv(pulse,x_upsample);
        
    case 2 %RRC
        span = 6;
        [pulse, t] = rtrcpuls(0.8,1/fsymb,fs,span);
        s = conv(pulse,x_upsample);
        
    otherwise
        disp('Unknown pulse');
        return;
end
t_vec = 1*Tsamp*(0:1:length(s)-1);
figure; subplot(2,1,1); plot(t_vec,real(s), 'b');
title('real')
subplot(2,1,2); plot(t_vec,imag(s), 'b');
title('imag')

% END OF EXERCISE 5
%% EX. 6

% Add AWGN noise to the signal
SNRdB = 10; %decide noise level; try snr=5
y = awgn(s, SNRdB, 'measured'); % add noise


MF = fliplr(conj(pulse));        %create matched filter impulse response
% figure;plot(MF);hold on; plot(pulse)
MF_output = conv(MF, y);  % Run through matched filter

figure(3); clf
subplot(2,1,1)
plot(real(s),'Linewidth',1); hold on
plot(real(y(1:end)),'Linewidth',1);
plot(real(MF_output),'Linewidth',1);
plot(real(MF_output(span*fsfd:end-span*fsfd)),'Linewidth',2)
legend('Tx signal', 'AWGN output', 'MF output with transients', 'MF output w/o transients')

% Rectangular pulse
if pulse_flag==0
    MF_output_cut = MF_output(2*span*fsfd : end-2*span*fsfd+1);%remove transient in the beginning and the end
else
    % RC or RRC pulse
    MF_output_cut = MF_output(2*span*fsfd-1  : end-2*span*fsfd+1);%remove transient in the beginning and the end
end
rx_vec = MF_output_cut(1:fsfd:end);  % downsampling

scatterplot(rx_vec); %scatterplot of received symbols

figure(3)
subplot(2,1,2);
plot(real(MF_output_cut));
hold on; stem(fsfd*(0:length(rx_vec)-1)+1,real(rx_vec));

% Minimum Eucledian distance detector
% Relate the detection to Detection region
metric = abs(repmat(rx_vec.',1,4) - repmat(const, length(rx_vec), 1)).^2; % compute the distance to each possible symbol
[tmp, m_hat] = min(metric, [], 2); % find the closest for each received symbol
m_hat = m_hat'-1;   % get the index of the symbol in the constellation

SER = sum(m-1 ~= m_hat) %count symbol errors
m_hat = de2bi(m_hat, 2, 'left-msb')'; %make symbols into bits
a_hat = m_hat(:)'; %write as a vector
BER = sum(a ~= a_hat) %count of bit errors


%% EX. 7 
% with time synchronization error
% eyediagram(MF_output, fsfd); % plot the eyediagram from the output of matched filter using MATLAB's function
% without time synchronization error
eyediagram(MF_output_cut, fsfd); % plot the eyediagram from the output of matched filter using MATLAB's function

% 
% %HOW to make your own eyediagram
%MF_eye = reshape(MF_output_cut(1:end), fsfd, []); %reshape in terms of samples, note that we will be fsfd/2 samples shifted
%MF_eye = circshift(MF_eye, fsfd/2);
%figure; plot(real(MF_eye), 'b'); %plot eye diagram, shift fsfd/2 to find sampling instant
% 
% 
% % BONUS: What happens to the eye diagram when there is a phase shift?
% % Here, we see that if we still want to find the correct sampling instant,
% % just take the modulus and we get rid of the phase.
phase_offset = exp(pi/3j);
MF_output_cut= MF_output_cut.*phase_offset;
eyediagram(MF_output_cut, fsfd); %phase error
% scatterplot(rx_vec*phase_offset)
% eyediagram(abs(MF_output_cut), fsfd);



