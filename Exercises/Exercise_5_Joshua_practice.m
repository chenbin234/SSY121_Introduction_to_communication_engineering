clc;close all

f_sample = 44000;
T_sample = 1/f_sample;
Rb = 440;
N = 432; % number of bits to tramsmit
fc = 1000; % carrier frequency

const = [(1+1i),(1-1i),(-1-1i),(-1+1i)]/sqrt(2);
% scatterplot(const);
% grid on;

M = length(const);
bits_per_symb = log2(M); % number of bits per symbol

Rs = Rb/bits_per_symb;
Ts = 1/Rs;
number_of_sample_per_symb = f_sample/Rs; % fsfd

% 随机生成一个消息
a = randsrc(1,N,[0,1]);
m = buffer(a,bits_per_symb)';
m_idx = bi2de(m, 'left-msb')+1;  % change m from binary to decimal(+1)
x = const(m_idx);  % 将m里面的每个symbol和const里面的点对应起来

x_upsample = upsample(x,number_of_sample_per_symb);

% Rectangular pulse
% pulse = ones(1,number_of_sample_per_symb);
% pulse = pulse./norm(pulse);

% RRC
span = 6;
[pulse, ~]= rtrcpuls(0.8,1/Rs,f_sample,span);
figure;
plot(pulse);
hold on;
%scatter(pulse);
title('pulse的时域图像')
length_of_pulse = length(pulse);

s = conv(pulse,x_upsample);


figure;
subplot(2,1,1)
t_vec = 1*T_sample*(0:1:length(s)-1);
plot(t_vec,real(s),'b');

samples = s(span*number_of_sample_per_symb:number_of_sample_per_symb:end-span*number_of_sample_per_symb);
hold on;
t = span*number_of_sample_per_symb:number_of_sample_per_symb:(span*number_of_sample_per_symb+number_of_sample_per_symb*(length(samples)-1));
hold on;
stem(t*T_sample,real(samples),'r')
title('real')
xlabel('seconds')

subplot(2,1,2);
plot(t_vec,imag(s),'b');
title('imag')
xlabel('seconds')

tx_signal = s.*exp(-1i*2*pi*fc*(0:length(s)-1)*T_sample);
tx_signal = real(tx_signal);
tx_signal = tx_signal/max(abs(tx_signal)); 

figure;
plot(tx_signal);


%% receiver part

% Add AWGN noise to the signal
SNRdB = 5; %decide noise level; try snr=5
%y = awgn(s, SNRdB, 'measured'); % add noise
y = awgn(tx_signal, SNRdB, 'measured'); % add noise
y = y.*exp(1i*2*pi*fc*(0:length(s)-1)*T_sample);

MF = fliplr(conj(pulse)); %create matched filter impulse response
MF_output = conv(MF, y);  % Run through matched filter

figure(3);
subplot(2,1,1)
plot(real(tx_signal),'Linewidth',1); hold on
plot(real(y(1:end)),'Linewidth',1);
plot(real(MF_output),'Linewidth',1);
plot(real(MF_output(span*number_of_sample_per_symb:end-span*number_of_sample_per_symb)),'Linewidth',2)
legend('Tx signal', 'AWGN output', 'MF output with transients', 'MF output w/o transients')

MF_output_cut = MF_output(2*span*number_of_sample_per_symb-1  : end-2*span*number_of_sample_per_symb+1); %卷积后的信号前面部分和后面部分不需要，需要剪掉
rx_vec = MF_output_cut(1:number_of_sample_per_symb:end); % 清除前面在x里面插入的0
scatterplot(rx_vec);

figure(3)
subplot(2,1,2);
plot(real(MF_output_cut));
hold on; 
stem(number_of_sample_per_symb*(0:length(rx_vec)-1)+1,real(rx_vec));

% Minimum Eucledian distance detector
% Relate the detection to Detection region
metric = abs(repmat(rx_vec.',1,4) - repmat(const, length(rx_vec), 1)).^2; % compute the distance to each possible symbol
[tmp, m_hat] = min(metric, [], 2); % find the closest for each received symbol
m_hat = m_hat'-1;   % get the index of the symbol in the constellation


SER = sum(m_idx'-1 ~= m_hat) %count symbol errors
m_hat = de2bi(m_hat, 2, 'left-msb')'; %make symbols into bits
a_hat = m_hat(:)'; %write as a vector
BER = sum(a ~= a_hat) %count of bit errors

