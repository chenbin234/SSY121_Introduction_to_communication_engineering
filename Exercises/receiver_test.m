clc;close all

f_sample = 5000; %how many samples do you pick per second
T_sample = 1/f_sample; % how long it take to pick a sample
Rb = 250; % how many bits you send per second
N = 432; % number of bits to tramsmit
fc = 1000; % carrier frequency

const = [(1+1i),(1-1i),(-1+1i),(-1-1i)]/sqrt(2);% need to be Gray code, be careful with the position
M = length(const); % how many points in the constellation
bits_per_symb = log2(M); % number of bits per symbol

Rs = Rb/bits_per_symb; % how many symbols you send per second
Ts = 1/Rs; % the time interval between symbolA and the after symbolB
fsfd = f_sample/Rs; % the number of sample in 1 Ts long

% 随机生成一个消息
a = randsrc(1,N,[0,1]);
m = buffer(a,bits_per_symb)'; % 将a切割成一个一个symbol
m_idx = bi2de(m, 'left-msb')+1;  % change m from binary to decimal(+1)
x = const(m_idx);  % 将m里面的每个symbol和const里面的点对应起来

preamble = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 ]; % Preamble, 13 bits from Barker code
x = [preamble, x];

x_upsample = upsample(x,fsfd);


% RRC
span = 6;
[pulse, ~]= rtrcpuls(0.4,1/Rs,f_sample,span);
% [pulse, ~]= rtrcpuls(0.4,3*(1/Rs),f_sample,span);
length_of_pulse = length(pulse);

s = conv(pulse,x_upsample);
figure;
subplot(2,1,1)
t_vec = 1*T_sample*(0:1:length(s)-1);
plot(t_vec,real(s),'b');

samples = s(span*fsfd:fsfd:end-span*fsfd);
hold on;
t = span*fsfd:fsfd:(span*fsfd+fsfd*(length(samples)-1));
hold on;
stem(t*T_sample,real(samples),'r')
title('real')
xlabel('seconds')

subplot(2,1,2);
% plot(t_vec,imag(s),'b');
plot(t_vec,imag(s),'b');
title('imag')
xlabel('seconds')

tx_signal = s.*exp(-1i*2*pi*fc*(0:length(s)-1)*T_sample);
tx_signal = real(tx_signal);
tx_signal = tx_signal/max(abs(tx_signal)); 

figure;
plot(t_vec,tx_signal);


%% receiver part

% Add AWGN noise to the signal
SNRdB = 100; %decide noise level; try snr=5
y = awgn(tx_signal, SNRdB, 'measured'); % add noise
y = y.*exp(1i*2*pi*fc*(0:length(tx_signal)-1)*T_sample); % y对应transmitter里面的s

% y(48199) = 100+100i;
% y(48200) = 1000+100i;
% y(48201) = 100+100i;
% y(48202) = 1000-100i;

y = y./max(abs(y));

preamble_upsample = upsample(preamble, fsfd);         % upsample preamble
conv_preamble_pulse = conv(pulse, preamble_upsample); % pulse shaping the upsampled preamble
corr = conv(y, fliplr(conv_preamble_pulse));  % correlate the recorded signal(maybe have the message) with processed preamble

corr = corr./13;                          % normalize corr 
plot(abs(corr))



sample_y = y(span*fsfd:fsfd:end-span*fsfd);

MF = fliplr(conj(pulse)); %create matched filter impulse response
MF_output = conv(MF, y);  % Run through matched filter


figure(3);
subplot(2,1,1)
plot(real(tx_signal),'Linewidth',1); hold on
plot(real(y(1:end)),'Linewidth',1);
plot(real(MF_output),'Linewidth',1);
plot(real(MF_output(span*fsfd:end-span*fsfd)),'Linewidth',2)
legend('Tx signal', 'AWGN output', 'MF output with transients', 'MF output w/o transients')

MF_output_cut = MF_output(2*span*fsfd-1  : end-2*span*fsfd+1); %卷积后的信号前面部分和后面部分不需要，需要剪掉
rx_vec = MF_output_cut(1:fsfd:end); % 清除前面在x里面插入的0

rx_vec = rx_vec(1+length(preamble):end);
rx_vec = rx_vec./abs(rx_vec);
scatterplot(rx_vec);

figure(3)
subplot(2,1,2);
plot(real(MF_output_cut));
hold on; 
stem(fsfd*(0:length(rx_vec)-1)+1,real(rx_vec));

% Minimum Eucledian distance detector
% Relate the detection to Detection region
metric = abs(repmat(rx_vec.',1,4) - repmat(const, length(rx_vec), 1)).^2; % compute the distance to each possible symbol
[tmp, m_hat] = min(metric, [], 2); % find the closest for each received symbol
m_hat = m_hat'-1;   % get the index of the symbol in the constellation


SER = sum(m_idx'-1 ~= m_hat) %count symbol errors
m_hat = de2bi(m_hat, 2, 'left-msb')'; %make symbols into bits
a_hat = m_hat(:)'; %write as a vector
BER = sum(a ~= a_hat) %count of bit errors

