% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency

function transmitter(pack,fc)

% Input parameters -–––––––––––––––––––––––––––––––––––––––––––––––––––––––
fs = 10000;           % Sampling frequency [Hz] (necessary for passband implementation)
Tsamp = 1/fs;         % Sampling time
Rb = 250;             % Bit rate [bit/sec] %Rb = fsymb*bpsymb; % Bit rate [bit/s]
N = 432;              % Number of bits to transmit
% fc = 1000;          % Carrier frequency [Hz]

const = [(1+1i),(1-1i),(-1+1i),(-1-1i)]/sqrt(2); % Constellation: QPSK/4-QAM

M = length(const);    % Number of symbols in the constellation
bpsymb = log2(M);     % Number of bits per symbol
fsymb = Rb/bpsymb;    % Symbol rate [symb/s]
Tsymb = 1/fsymb;      % Symbol time
fsfd = fs/fsymb;      % Number of samples per symbol [samples/symb] (choose fs such that fsfd is an integer for simplicity)

alpha = 0.8;          % Roll off factor / Excess bandwidth factor (a_RC=0.35;a_RRC=0.8)
tau = 1/fsymb;        % Nyquist period or symbol time 
span = 6;             % Pulse width (symbol times of pulse)

% Bit to symbol mapping & spacing: -–––––––––––––––––––––––––––––––––––––––

m = buffer(pack, bpsymb)';                    % Group bits into bits per symbol
m_idx = bi2de(m, 'left-msb')'+1;              % Bits to symbol index
x = const(m_idx);                             % Look up symbols using the indices

% Add preamble: -––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
preamb = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];     % 13 bits from Barker code
x = [preamb, x]; 

x_upsample = upsample(x, fsfd);               % Space the symbols fsfd apart, to enable pulse shaping using conv.

% Pulse shaping - Convert symbols to baseband signal: -––––––––––––––––––––

[pulse,~] = rtrcpuls(alpha,tau,fs,span);      % Create rrc pulse: rtrcpuls(alpha,tau,fs,span)
s = conv(pulse,x_upsample);                   % Create baseband signal (convolve symbol with pulse)


% Add carrier - Convert into passband signal: -––––––––––––––––––––––––––––

tx_signal = s.*exp(-1i*2*pi*fc*(0:length(s)-1)*Tsamp); % Put the baseband signal on a carrier signal

% Modulation/Upconversion: -–––––––––––––––––––––––––––––––––––––––––––––––
tx_signal = real(tx_signal);                  % Send real part, information is in amplitude and phase
tx_signal = tx_signal/max(abs(tx_signal));    % Limit the max amplitude to 1 to prevent clipping of waveforms
disp(['the length of transmitter signal is',num2str(length(tx_signal))])


% Transmission: -––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

player = audioplayer(tx_signal, fs); % Create an audioplayer object to play the noise at a given sampling frequency
playblocking(player); % Play the noise 