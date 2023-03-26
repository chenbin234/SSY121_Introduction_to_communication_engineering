func3 =@(alpha,tau,tvec) sinc(tvec/tau).*cos(alpha*pi*tvec/tau)./(1-(2*alpha*tvec/tau).^2);
% func3 =@(alpha,tau,tvec) sinc(tvec/tau).*cos((pi*alpha/tau)*tvec)./(1-((2*alpha/tau)*tvec).^2);
T_symb = 0.01;
alpha = 0.3;
BW = (1+alpha)/(2*T_symb);

f_samp = 2*BW;
span = 6;

tvec = 0:(1/f_samp):span*T_symb;
tvec = [-fliplr(tvec(2:end)) tvec];

fss = 50*BW;
cont_tvec = eps:(1/fss):span*T_symb;
cont_tvec = [-fliplr(cont_tvec(2:end)) cont_tvec];

p_samp = func3(alpha,T_symb,tvec);
p_continuous = func3(alpha,T_symb,cont_tvec);

figure;
plot(cont_tvec,p_continus,'r','LineWidth',2);
hold on;
stem(tvec,p_samp);
xlabel('t(s)');
legend('Continus signal','Sampled signal');
title(['Sampled signal, f_{samp}/BW= ',num2str(f_samp/BW)]);

%% Perform DFT and convert the digital frequence to corresponding analog
N = length(p_samp);
P = fftshift(fft(p_samp));
phase = unwrap(angle(P));
df = f_samp/N;
fvec = df*(-floor(N/2):1:ceil(N/2)-1);

figure;
subplot(2,1,1);
plot(fvec,abs(P)/N,'-o','LineWidth',2);
title(['Frequency response of sampled signal, f_{samp}/BW= ',num2str(f_samp/BW)]);
xlabel('Frequency in Hz');
ylabel('Power in dB');
xlim([-1*BW,1*BW]);

subplot(2,1,2);
plot(fvec,phase*180/pi,'LineWidth',2);
title(['Phase response of sampled signal, f_{samp}/BW= ',num2str(f_samp/BW)]);
xlabel('Frequency in Hz')
ylabel('Phase (Degrees)')
xlim([-1*BW,1*BW])

%% reconstruct
[p_recons,t_recons] = reconstruct_sinc(p_samp,f_samp);
figure;
plot(t_recons,p_recons,'LineWidth',2);
hold on;
plot(cont_tvec,p_continuous,'r','LineWidth',2);
title(['Reconstructed signal, f_{samp}/BW= ',num2str(f_samp/BW)]);
legend('Reconstructed','Equivalent to Continuous')

function [p_recons,t_recons] = reconstruct_sinc(p_samp,f_samp)

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

















