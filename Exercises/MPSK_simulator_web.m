
% System settings
Q=9;   % # samples per symbol, Q=(sampling rate)/(symbol rate)
P=-13; % floor(abs(P)) is the # symbols in the PA, negative P creates a Barker PA of length {2,3,4,5,7,11,13}, integer P>0 creates a BPSK PA, non-integer P>1.0 creates Gaussian PA, imaginary P create complex-valued PA
D=10;  % # data symbols
M=2;   % # of symbol alternatives in M-PSK

% Channel settings
delay=randi(50,1)*1; % random delay in number of samples
delay_offset=0;      % delay offset in number of samples
phi=randi(360,1)*1;  % random phase between 0 and 360 degrees
phi_offset=0;        % phase offset in degrees
freq_offset=0;       % degrees per symbol
sigma=0.2;           % standard deviation of the AWGN

% Basic pulse
v=ones(1,Q);     % rectangular
%v=ones(1,Q+2);  % rectangular with ISI
%v=(Q-1)/2-abs(-(Q-1)/2:(Q-1)/2); % triangle
%v=1:Q;          % ramp
%v=randn(1,Q);   % random

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FontSize=13;
for k=1:3, figure(k), clf, end

% Transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

% Create Preamble (PA)
if (real(P)<0)||(imag(P)<0), switch round(abs(P)), case 2, PA_symb=[+1 -1]; case 3, PA_symb=[+1 +1 -1]; case 4, PA_symb=[+1 +1 -1 +1]; case 5, PA_symb=[+1 +1 +1 -1 +1]; case 7, PA_symb=[+1 +1 +1 -1 -1 +1 -1]; case 11, PA_symb=[+1 +1 +1 -1 -1 -1 +1 -1 -1 +1 -1]; case 13, PA_symb=[+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1]; otherwise, PA_symb=sign(randn(1,P)); end, if ~isreal(P), PA_symb=(PA_symb-j*PA_symb)/sqrt(2); end, elseif mod(abs(P),1)==0, if isreal(P), PA_symb=sign(randn(1,abs(P))); else, PA_symb=(sign(randn(1,abs(P)))+j*sign(randn(1,abs(P))))/sqrt(2); end, else, if isreal(P), PA_symb=randn(1,round(abs(P))); else, PA_symb=randn(1,round(abs(P)))+j*randn(1,round(abs(P))); end, PA_symb=PA_symb/sqrt(sum(abs(PA_symb).^2)/round(abs(P))); end, P=floor(abs(P));

% Data
d=randi(M,1,D)-1; % random data symbol index {0:(M-1)}
D_symb=exp(j*2*pi*d/M+j*pi/4*(M==4)); % M-PSK constellation
PAD_symb=[PA_symb, D_symb, zeros(1,P+D)];
subplot(411), stem(real(PAD_symb),'b','linewidth',2), hold on, if ~isreal(PA_symb)||(M>2), stem(imag(PAD_symb),'r','linewidth',2), end, plot([1 1]*(P+0.5),[-1 1]*2,'k--','linewidth',2), hold off, if ~isreal(PA_symb)||(M>2), legend('real','imag','Location','NorthEast'), end, ylabel('PA & Data symbols'), axis([0 P+D+1 [-1 1]*max(abs(PAD_symb))*1.1]), title(sprintf('Q = %d, P = %d, D = %d, M = %d',Q,P,D,M)), grid on, set(gca,'FontSize',FontSize)

% Train
dirac=[1 zeros(1,Q-1)];
PAD_train=kron(PAD_symb,dirac);
subplot(412), stem(real(PAD_train),'b','linewidth',2), hold on, if ~isreal(PA_symb)||(M>2), stem(imag(PAD_train),'r','linewidth',2), end, plot([1 1]*(P*Q-0.5*Q+1),[-1 1]*2,'k--','linewidth',2), hold off, if ~isreal(PA_symb)||(M>2), legend('real','imag','Location','NorthEast'), end, ylabel('PA & Data train'), axis([0 (P+D)*Q+1 [-1 1]*max(abs(PAD_train))*1.1]), grid on, set(gca,'FontSize',FontSize)

% Basic pulse
v=v/sqrt(sum(v.^2)); % normalize to energy = 1
subplot(413), plot(1:length(v),v,'b',[1 1 Q Q],[-1 1 1 -1]*2,'k--','linewidth',2), ylabel('Basic pulse v(t)'), axis([0 length(v)+1 [min([v,-0.01]) max(v)]*1.1]), grid on, set(gca,'FontSize',FontSize)

% TX signal
s=conv(PAD_train,v); % TX filter
subplot(414), plot(real(s),'b','linewidth',2), hold on, if ~isreal(PA_symb)||(M>2), plot(imag(s),'r','linewidth',2), end, plot([1 1]*(P*Q+0.5),[-1 1],'k--','linewidth',2), hold off, if ~isreal(PA_symb)||(M>2), legend('real','imag','Location','NorthEast'), end, ylabel('TX signal s(t)'), grid on, axis([0 (P+D+1)*Q [-1 1]*max(abs(s))*1.1]), set(gca,'FontSize',FontSize)

% Channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

% Frequency shift
s=s.*exp(j*freq_offset/180*pi/Q*[[0:length(s)-1]-P*Q/2+0.5]);

% Phase shift
s=s*exp(j*phi/180*pi);

% Delay
s=[zeros(1,delay), s];

% AWGN
y=s+(randn(size(s))+j*randn(size(s)))*sigma;
subplot(411), plot(1:length(y),real(y),'b',1:length(y),imag(y),'r','linewidth',2), hold on, plot([1 1]*(P*Q+delay+0.5),[-1 1]*9,'k--','linewidth',2), hold off, ylabel('RX signal y(t)'), legend('real','imag','Location','NorthEast'),  grid on, axis([0 (P+D+1)*Q+delay [-1 1]*max(abs(y))*1.1]), set(gca,'FontSize',FontSize), title(sprintf('\\Delta = %d samples, \\phi = %2.1f\\circ, \\delta = %2.1f\\circ/symbol, \\sigma = %2.1f',delay,phi,freq_offset,sigma)), set(gca,'FontSize',FontSize)

% Receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PA correlation
PA_train=kron(PA_symb,dirac);
s_PA=conv(PA_train,v);
z_PA=conv(y,fliplr(conj(s_PA)));

% Detect the PA and estimate delay and phase
[peak, ind]=max(abs(z_PA));
ind=ind+delay_offset;
if P>0 % if a PA is transmitted
    delay_hat=ind+1-Q*(P+1);
    phi_hat=mod(angle(z_PA(ind))*180/pi+phi_offset,360);
    subplot(412), plot(1:length(z_PA),real(z_PA),'b',1:length(z_PA),imag(z_PA),'r','linewidth',2), ylabel('PA correlation'), legend('real','imag','Location','NorthEast'), grid on, axis([0 (2*P+D+2)*Q+delay-2 [-1 1]*max(0.01+abs([real(z_PA), imag(z_PA)]))*1.1]), set(gca,'FontSize',FontSize)
    subplot(413), plot(1:length(z_PA),abs(z_PA),'b',ind,abs(z_PA(ind)),'r*','linewidth',2), ylabel('abs(PA correlation)'), grid on, axis([0 (2*P+D+2)*Q+delay-2 [0 1]*max(0.01+abs(z_PA))*1.1]), title(sprintf('\\Delta = %d samples, \\Delta* = %d samples',delay,delay_hat)), set(gca,'FontSize',FontSize), set(gca,'YTick',unique([0 round(max(abs(z_PA((ind-Q):-Q:1))),1) P]))
else % if no PA is included
    delay_hat=delay+delay_offset;
    phi_hat=mod(phi+phi_offset,360);
end

% Matched filter
data_start=delay_hat+Q*(P+1);
data_ind=data_start:Q:(data_start+(D-1)*Q);
z=conv(y,fliplr(v));
subplot(414), plot(1:length(z),real(z),'b',1:length(z),imag(z),'r',data_ind,real(z(data_ind)),'b*',data_ind,imag(z(data_ind)),'r*','linewidth',2), hold on, plot([1 1]*(P*Q+delay+0.5*Q),[-1 1]*9,'k--','linewidth',2), hold off, ylabel('MF output z(t)'), legend('real','imag','Location','NorthEast'), grid on, axis([0 (1*P+D+2)*Q+delay-1 [-1 1]*max(abs([real(z), imag(z)]))*1.1]), set(gca,'FontSize',FontSize)

% Diagnostics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

% Phase compensation
z_rot=z*exp(-j*phi_hat/180*pi);

% Decisions on the symbols and error rates
d_hat=mod(round((angle(z_rot(data_ind))-pi/4*(M==4))/(pi*2)*M),M);
Symb_errors=sum(d~=d_hat);
SER=Symb_errors/D;
SER_UB=(1+(M>2))*qfunc(sin(pi/M)/sigma); % Upper bound

if (P>0)&&(D>0)
    subplot(221), plot([0 real(exp(j*phi/180*pi))]*abs(peak)*2,[0 imag(exp(j*phi/180*pi))]*abs(peak)*2,'k--',real(z_PA),imag(z_PA),'b',real(z_PA(ind)),imag(z_PA(ind)),'r*','linewidth',2), axis square, axis(abs(peak)*1.5*[-1 1 -1 1]), grid on, xlabel('real'), ylabel('imag'), title(sprintf('PA correlation, \\phi = %2.1f\\circ, \\phi* = %2.1f\\circ',phi,phi_hat)), set(gca,'FontSize',FontSize)
elseif P>0
    subplot(111), plot([0 real(exp(j*phi/180*pi))]*abs(peak)*2,[0 imag(exp(j*phi/180*pi))]*abs(peak)*2,'k--',real(z_PA),imag(z_PA),'b',real(z_PA(ind)),imag(z_PA(ind)),'r*','linewidth',2), axis square, axis(abs(peak)*1.5*[-1 1 -1 1]), grid on, xlabel('real'), ylabel('imag'), title(sprintf('PA correlation, \\phi = %2.1f\\circ, \\phi* = %2.1f\\circ',phi,phi_hat)), set(gca,'FontSize',FontSize)
end

if D>0
    z_eye=reshape(z_rot((data_start-floor(Q/2)):(data_start-floor(Q/2)+D*Q-1)),Q,D);
    subplot(222), plot((0:Q-1)-floor(Q/2),real(z_eye(:,1)),'b','linewidth',2), hold on, if M<=2, plot((0:Q-1)-floor(Q/2),real(z_eye),'b','linewidth',2), else, plot((0:Q-1)-floor(Q/2),imag(z_eye),'r','linewidth',2), plot((0:Q-1)-floor(Q/2),real(z_eye),'b','linewidth',2), end, plot([0 0],[-1 1]*9,'k--','linewidth',2), if M>2, legend('real','imag','Location','NorthEast'), end, hold off, axis square, axis([-floor(Q/2) Q-1-floor(Q/2) [-1 1]*max(abs(z_rot))*1.1]), grid on, xlabel('sample index'), title(sprintf('Eye diagram, \\Delta = %d samples, \\Delta* = %d samples',delay,delay_hat)), set(gca,'FontSize',FontSize)
    
    subplot(223), circ=exp(j*(0:0.01:2*pi))*max(sigma,0.03); ColorOrder=[0.00 0.00 1.00; 1.00 0.00 0.00; 0.13 0.94 0.96; 0.11 0.78 0.39; 0.80 0.43 0.91; 0.97 0.76 0.29; 0.74 0.27 0.25; 0.20 0.31 0.20]; for k=0:M-1, plot(real(z(data_ind(d==k))),imag(z(data_ind(d==k))),'*','Color',ColorOrder(mod(k,8)+1,:),'linewidth',2), hold on, end, for k=0:M-1, region=[0 exp(j*2*pi*(k+0.5)/M+j*pi/4*(M==4))]*9; plot(real(region),imag(region),'k--','linewidth',2), plot(real(circ+exp(j*2*pi*k/M+j*pi/4*(M==4))),imag(circ+exp(j*2*pi*k/M+j*pi/4*(M==4))),'Color',ColorOrder(mod(k,8)+1,:),'linewidth',2), end, hold off, axis square, axis((1+max(sigma,0.1)*3)*[-1 1 -1 1]), grid on, xlabel('real'), ylabel('imag'), title('MF output sampled per symbol'), set(gca,'FontSize',FontSize)
    
    subplot(224), circ=exp(j*(0:0.01:2*pi))*max(sigma,0.03); ColorOrder=[0.00 0.00 1.00; 1.00 0.00 0.00; 0.13 0.94 0.96; 0.11 0.78 0.39; 0.80 0.43 0.91; 0.97 0.76 0.29; 0.74 0.27 0.25; 0.20 0.31 0.20]; for k=0:M-1, plot(real(z_rot(data_ind(d==k))),imag(z_rot(data_ind(d==k))),'*','Color',ColorOrder(mod(k,8)+1,:),'linewidth',2), hold on, end, for k=0:M-1, region=[0 exp(j*2*pi*(k+0.5)/M+j*pi/4*(M==4))]*9; plot(real(region),imag(region),'k--','linewidth',2), plot(real(circ+exp(j*2*pi*k/M+j*pi/4*(M==4))),imag(circ+exp(j*2*pi*k/M+j*pi/4*(M==4))),'Color',ColorOrder(mod(k,8)+1,:),'linewidth',2), end, hold off, axis square, axis((1+max(sigma,0.1)*3)*[-1 1 -1 1]), grid on, xlabel('real'), ylabel('imag'), title('De-rotated MF output sampled per symbol'), text(-1,1.1,sprintf('SER = %d / %d = %3.1f%% (%3.1f%%)',Symb_errors,D,SER*100,SER_UB*100),'FontSize',FontSize,'BackgroundColor','w'), set(gca,'FontSize',FontSize)
end
