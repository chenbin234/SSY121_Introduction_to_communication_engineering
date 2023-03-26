% Bandlimited signals made of sinc
func1 = @(Ts,tvec) 6*sinc(tvec/(2*Ts)).^2 + sinc(tvec/(6*Ts)).^6 + 6*sinc(tvec/(10*Ts)).^10; 
%% question(a)
T_symbol = 0.025;
BW = 1/(2*T_symbol);
span = 6;
f_sample = 2*BW;
% 采样时间间隔
tvec = 0:(1/f_sample):span*T_symbol;
tvec = [-fliplr(tvec(2:end)) tvec];
% 连续时间间隔
fss = 50*BW;
cont_tvec = eps:(1/fss):span*T_symbol;
cont_tvec = [-fliplr(cont_tvec(2:end)) cont_tvec];
% 调用函数
p_sample = func1(T_symbol,tvec);
p_continuous = func1(T_symbol,cont_tvec);
% 画图
figure;
plot(cont_tvec,p_continuous,'r','LineWidth',2);
hold on;
stem(tvec,p_sample);
xlabel('t(s)')
legend('Continuous signal','Sampled signal')
title(['Sampled signal, f_{samp}/BW= ', num2str(f_sample/BW)])

%% question(b)
N = length(p_sample);
P = fftshift(fft(p_sample));
phase = unwrap(angle(P));
df = f_sample/N;
fvec = df*(-floor(N/2):1:ceil(N/2)-1);

figure;
% 图1
subplot(2,1,1);
plot(fvec, abs(P)/N,'o-','LineWidth',2)
title(['Frequence response of sampled signal, f_{samp}/BW]= ',num2str(f_sample/BW)])
xlabel('Frequence in Hz')
ylabel('Power in dB')
xlim([-1*BW,1*BW])
% 图2
subplot(2,1,2);
plot(fvec,phase*180/pi,'LineWidth',2)
title(['Phase response of sampled signal, f_{samp}/BW= ',num2str(f_sample/BW)])
xlabel('Phase (Degrees)')
xlim([-1*BW,1*BW])

%% question(c)

%调用重建函数
[p_recons,t_recons] = reconstruct_sinc(p_sample,f_sample);

figure;
plot(t_recons,p_recons,'LineWidth',2);
hold on;
plot(cont_tvec,p_continuous,'r','LineWidth',2);
title(['Reconstructed signal, f_{samp}/BW= ',num2str(f_sample/BW)]);
legend('Reconstructed','Equivalent to Continuous')


function [p_recons,t_recons]=reconstruct_sinc(p_samp,f_samp)
span = 6;
cont_f_samp = 10*f_samp;
cont_t_samp = 1/cont_f_samp;

t_samp = 1/f_samp;
tvec = eps:cont_t_samp:span*t_samp;
tvec = [-fliplr(tvec(2:end)) tvec];
sinc_pulse = sinc(f_samp*tvec);

p_upsamp = upsample(p_samp,cont_f_samp/f_samp);
p_upsamp(end-cont_f_samp/f_samp+2:end)=[];
p_recons = conv(p_upsamp,sinc_pulse,'same');
t_recons = (1/cont_f_samp).*(-floor(length(p_recons)/2):1:floor(length(p_recons)/2));

end




