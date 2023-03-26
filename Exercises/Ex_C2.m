% Exercise C2

%% 5a)

clc, clear all, close all

m_size = 2; 
number_of_symbols = 10;
f_samp = 44000; % Sampling frequency
t_samp = 1/f_samp;
f_c = 1000;

Rb = 440;
Rs = Rb/log2(2^m_size);
Ts = 1/Rs;
sps = f_samp/Rs; % samples per symbol

a = randi([0 1], 1, number_of_symbols*m_size);

% Bits to messeges => Messeges to symbols 

symb = zeros(1, length(a)/m_size); % prep symb array
const = qammod(0:3,m_size^2);

% % Plot
% figure(1)
% scatter(real(const), imag(const), 'red', 'filled'), grid on
% axis([-3 3 -3 3])
% xlabel('Real axis')
% ylabel('Imaginary axis')

% Mapping of meeseges to symbols
for i = 1:length(symb)
    if (a(2*i-1)==0 && a(2*i)==0)
        symb(i)=const(3);
    elseif (a(2*i-1)==0 && a(2*i)==1)
        symb(i)=const(1);
    elseif (a(2*i-1)==1 && a(2*i)==0)
        symb(i)=const(4);
    elseif (a(2*i-1)==1 && a(2*i)==1)
        symb(i)=const(2);
    else
        print('Errrrrroooorr')
    end
end

figure(2)
scatter(real(symb), imag(symb), 'filled'), grid on
axis([-3 3 -3 3])
xlabel('Real axis')
ylabel('Imaginary axis')

% Symbols to waveform

symb_up = upsample(symb, sps);  % upsampling of signal

% Create pulse
[g, t] = rcpuls(0.5, t_samp, f_samp, 6);     % Create Raised Cosine
plot(t, g)

% Convolve with pulse
signal = conv(symb_up,g);

t1= (0:numel(signal)-1)*Ts;

plot(t1,signal), grid on


%% 










