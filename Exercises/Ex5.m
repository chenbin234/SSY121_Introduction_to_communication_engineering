% Signal space and constellation
% Ex 5
% By: Johan Oestman
clc; close all

% Input parameters
% For passband implemetation sampling frequency is a necessary parameter
fs = 44000; % sampling frequency [Hz]
Tsamp = 1/fs;    % Sampling time
Rb = 440; % bit rate [bit/sec]
N = 100; % number of bits to transmit
fc = 1000;
%----------------------------------------------
%       SET THIS FLAG TO CHOOSE WHAT PULSE! 
%           0 = rect, 1 = rc, 2 = rrc
pulse_flag = 2; 
%----------------------------------------------


% Constellation or bit to symbol mapping
const = [(1 + 1i), (1 - 1i), (-1 -1i), (-1 + 1i)]/sqrt(2);% Constellation 1 - QPSK/4-QAM

scatterplot(const); grid on;                          % Constellation visualization

M = length(const);                                    % Number of symbols in the constellation
bpsymb = log2(M);                                     % Number of bits per symbol
fsymb = Rb/bpsymb;                                    % Symbol rate [symb/s]
Tsymb = 1/fsymb;                                      % symbol time
fsfd = fs/fsymb;                                      % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity) [samples/symb]
a = randsrc(1,N,[0 1]);                               % Information bits
m = buffer(a, bpsymb)';                               % Group bits into bits per symbol
m_idx = bi2de(m, 'left-msb')'+1;                      % Bits to symbol index
x = const(m_idx);                                     % Look up symbols using the indices

x_upsample = upsample(x, fsfd);                       % Space the symbols fsfd apart, to enable pulse shaping using conv.



% Now we have a signal ready to be shaped.
switch pulse_flag
    case 0 %Rectangular
        pulse = ones(1, fsfd);                   % Generate a rectangular pulse
        pulse = pulse ./ norm(pulse);
        s = conv(pulse,x_upsample);    % Each symbol replaced by the pulse shape and added
        
    case 1 %RC
        span = 6; %how many symbol times to we want of pulse (note that RC and RRC is not zero outside of symbol time!)
        [pulse, t] = rcpuls(0.35,1/fsymb,fs,span);
        s = conv(pulse,x_upsample);
        
    case 2 %RRC
        span = 6;
        [pulse, t] = rtrcpuls(0.4,1/fsymb,fs,span);
        s = conv(pulse,x_upsample);
        
    otherwise
        disp('Unknown pulse');
        return;
end

figure; 
subplot(2,1,1); 
plot(Tsamp*(0:(length(s)-1)), real(s), 'b');

samples = s(span*fsfd:fsfd:end-span*fsfd);
hold on;
t = span*fsfd:fsfd:(span*fsfd + fsfd*(length(samples)-1));
hold on;
stem(t*Tsamp, real(samples),'r')


title('real')
xlabel('seconds')
subplot(2,1,2); 
plot(Tsamp*(0:(length(s)-1)), imag(s), 'b');
title('imag')
xlabel('seconds')


if pulse_flag > 0
    tx_signal = s.*exp(-1i*2*pi*fc*(0:length(s)-1)*Tsamp); % Carrier Modulation/Upconversion 
    tx_signal = real(tx_signal);                                        % send real part, information is in amplitude and phase
    tx_signal = tx_signal/max(abs(tx_signal));          % Limit the max amplitude to 1 to prevent clipping of waveforms

    figure
    plot(tx_signal)

end



