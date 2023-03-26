% RECEIVER 
%
% This is the receiver structure that you will have to complete.
% The function: receiver(fc) is a setup function for the receiver. Here,
% the audiorecorder object is initialized (see help audiorecorder or
% MATLAB's homepage for more information about the object).
% 
% The callback function audioTimerFcn() is a callback function that is
% triggered on a specified time interval (here it is determined in the
% setup function, by the variable time_value)
% 
% Your task is to extend this code to work in the project!
%%

function [audio_recorder] = receiver(fc)

fs = 10000;                                      % Sampling frequency
audio_recorder = audiorecorder(fs,24,1);         % Create the recorder

audio_recorder.UserData.fs = fs;                 % store fs information into the recObject
audio_recorder.UserData.fc  = fc;                % store fc information into the recObject 

% Attach callback function:
time_value = 0.2;                                % How often the function should be called, in seconds
audio_recorder.UserData.time_value = time_value; % store time_value information into the recObject
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); % attach a function that should be called every second, the function that is called is specified below.

flag = 0;
audio_recorder.UserData.flag = flag;             % store time_value information into the recObject
audio_recorder.UserData.counter = 1;             % count how many times before detect a new preamble

disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(['This time, the time_value=',num2str(time_value)])

% ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
audio_recorder.UserData.receive_complete = 0; % this is a flag that the while loop in the GUI will check
audio_recorder.UserData.pack  = []; %allocate for data package
audio_recorder.UserData.pwr_spect = []; %allocate for PSD
audio_recorder.UserData.const = []; %allocate for constellation
audio_recorder.UserData.eyed  = []; %allocate for eye diagram

record(audio_recorder); % Start recording
end

% CALLBACK FUNCTION
% This function will be called every [time_value] seconds, where time_value
% is specified above. Note that, as stated in the project MEMO, all the
% fields: pwr_spect, eyed, const and pack need to be assigned if you want
% to get outputs in the GUI.

% So, let's see an example of where we have a pulse train as in Computer
% exercise 2 and let the GUI plot it. Note that we will just create 432
% random bits and hence, the GUI will not be able to decode the message but
% only to display the figures.
% Parameters in the example:
% f_s = 22050 [samples / second]
% R_s = 350 [symbols / second]
% fsfd = f_s/R_s [samples / symbol]
% a = 432 bits
% M = 4 (using QPSK as in computer exercise)

function audioTimerFcn(recObj, event, handles)

disp('Callback triggered')

% ###### Basic parameter ######

fc = recObj.UserData.fc;                            % Carrier frequency
fs = recObj.UserData.fs;                            % Sampling frequency [Hz]
Tsamp = 1/fs;                                       % Sampling time
Rb = 250;                                           % Bit rate [bit/sec] %Rb = fsymb*bpsymb; % Bit rate [bit/s]
N = 432;                                            % Number of bits to transmit
const = [(1+1i), (1-1i), (-1+1i), (-1-1i)]/sqrt(2); % Constellation: QPSK/4-QAM
M = length(const);                                  % Number of symbols in the constellation
bpsymb = log2(M);                                   % Number of bits per symbol
fsymb = Rb/bpsymb;                                  % Symbol rate [symb/s]
Tsymb = 1/fsymb;                                    % Symbol time
fsfd = fs/fsymb;                                    % Number of samples per symbol [samples/symb] (choose fs such that fsfd is an integer for simplicity)
alpha = 0.8;                                        % Roll off factor / Excess bandwidth factor (a_RC=0.35;a_RRC=0.8)
span = 6;                                           % Pulse width (symbol times of pulse)
% Bw = (1+alpha)/(2*Tsymb);                         % Bandwidth
preamble = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];         % Preamble, 13 bits from Barker code
[pulse,~] = rtrcpuls(alpha,Tsymb,fs,span);          % RRC pulse
disp('Parameters are set')


% ###### Detect the value of time_value ######

time_value = recObj.UserData.time_value;            % time_value when no message are found
flag = recObj.UserData.flag;


if flag == 0                                        % 1.means No preamble found yet 2. or means we finished processing last received message

    % ###### Filter the received signal ######
    disp(['The value of flag now is: ',num2str(recObj.UserData.flag)])
    disp(['When start to detect preamble, the time value is:',num2str(time_value)])
    
    rx = getaudiodata(recObj)';                     % Get record data from recObj
    rx = rx(end-(time_value*fs):end);          
    disp('Signal is captured')
 
    y = rx.*exp(1i*2*pi*fc*(0:length(rx)-1)*Tsamp); % Get baseband signal
    y = y./max(abs(y));                             % normalise y
    
    disp('Signal is demodulated')
    
    % In order to find out the start of message, barker code was used as preamble
    % What we've done to preamble in transmitter should be done a second time
    % in order to take advantage of its low autocorrelation property
    
    preamble_upsample = upsample(preamble, fsfd);         % upsample preamble  
    conv_preamble_pulse = conv(pulse, preamble_upsample); % pulse shaping the upsampled preamble
    corr = conv(y,fliplr(conv_preamble_pulse));           % correlate the recorded signal(maybe have the message) with processed preamble
    disp('correlate the recorded signal with conv_preamble_pulse')

    corr = corr./13;                                      % normalize corr 
    % figure(2); clf;plot(abs(corr));
    
    Threshold = 0.8;                                      % if corr has a peak over 0.8, there are messages in transmitter
    [tmp,Tmax] = max(abs(corr));                          % value and location of the peak of corr


    % if corr has a peak over threshold, there are messages in rx, otherwise there is
    % no message in rx

    if (tmp < Threshold)
        
        recObj.UserData.counter = recObj.UserData.counter + 1;
        disp('No (preamble)message find!');               % When tmp < Threshold, that means no preamble found, that is no message found
        disp(recObj.UserData.counter)
        disp(' ')
        disp(' ')
        disp(' ')
        disp(' ')
        time_value = 0.2;
        set(recObj,'TimerPeriod',time_value);
    else                                                  % When tmp > Threshold, that means (preamble)message is found
    
        disp('find the preamble!')
        resume(recObj);
        time_value = 2;                                   % set time_value to 2s in order to get the whole message 
        set(recObj,'TimerPeriod',time_value); 
        recObj.UserData.time_value = 2;
        disp(['When we first found the preamble, the time value is changed to:',num2str(time_value)])
        
        recObj.UserData.flag = 1;                         % set flag to 1 in order to process the message in next round

    end

else

    % it means flag=1, we found preamble in last round
    % we now try to get message and decode it.

    disp(['The value of flag now is: ',num2str(recObj.UserData.flag)])
    disp(['the time value now is: ',num2str(time_value)])
    disp(['the time value now is: ',num2str(recObj.UserData.time_value)])
    rx = getaudiodata(recObj)';                           % Get record data from recObj
    rx = rx(end-(2.4)*fs:end);                            % get the message signal+noise
        
    y = rx.*exp(1i*2*pi*fc*(0:length(rx)-1)*Tsamp);       % Get baseband signal
    y = y./max(abs(y));                                   % normalise y
    
    % Lowpass filter:
    [a,b]=butter(2,0.05);
    y=filtfilt(a,b,y);
    
    disp(["the coarse signal's (longer than we signal we want) length is )",num2str(length(y))])
    disp('The message is demodulated')
    
    preamble_upsample = upsample(preamble, fsfd);         % upsample preamble  
    conv_preamble_pulse = conv(pulse, preamble_upsample); % pulse shaping the upsampled preamble
    corr = conv(y,fliplr(conv_preamble_pulse));           % correlate the recorded signal(maybe have the message) with processed preamble
    
    [tmp,Tmax] = max(abs(corr));                          % value and location of the peak of corr
    
    Tx_hat = Tmax - length(conv_preamble_pulse);          % find delay, Tx_hat is the location of the start(preamble) of the message
    display(['The Tmax is ',num2str(Tmax)])
    display(['The length of conv_preamble_pulse is ',num2str(length(conv_preamble_pulse))])
    display(['The Tx_hat is ',num2str(Tx_hat)])
 
    length_signal = (fsfd*(length(preamble)+N./bpsymb)+length(pulse)-1); %length of preamble+message in y
    disp(['The theoretical length of signal we should capture is',num2str(length_signal)])

    y = y((Tx_hat+1):(Tx_hat+length_signal));             % get the segment of y we want - which is preamble+message
    disp(['Try to capture the whole signal (preamble+message), the length we captured is:',num2str(length(y))])
    disp('the message is located')
    
    % Matched filter
    MF = fliplr(conj(pulse));                             % Create matched filter impulse response
    MF_output = conv(MF, y);                              % Run through matched filter
    
    disp('The signal is filtered')
    
    MF_output_cut = MF_output(2*span*fsfd-1:end-2*span*fsfd+1);  % cut the beginning and end of match filter output
    MF_output_cut_without_premable = MF_output(2*span*fsfd-1+length(preamble)*fsfd:end-2*span*fsfd+1); % cut the preamble part
    rx_preamble_message = MF_output_cut(1:fsfd:end);      % dowmsampling, get the preamble+message                 
    rx_vec = rx_preamble_message(1+length(preamble):end); % get the message
    
        
    % Phase synchronization & Frequency synchronization
    rx_preamble = rx_preamble_message(1:length(preamble));% extract the recieved preamble symbols
    rx_preamble = rx_preamble./mean(abs(rx_preamble));    % normalize the preamle symbols
    
    phase = (angle(rx_preamble) - angle(preamble)); % calculated the difference in angle between known preamble and received preamble  
    disp('Phase offset before mod:')   
    disp(phase)
    phase = mod(phase, 2*pi); % Get offset  
    
    % angle(rx_preamble) - angle(preamble) = phase_offset + 2*pi*frequency_offset 
    % polyfit(x,phase), the intercept is phase offset, the slope is 2*pi*frequency_offset
    % we can got the phase offset and frequency_offset
    x = Tsamp*[1:1:length(preamble)];
    P = polyfit(x,phase,1);
    slope = P(1);
    fdelta = slope./(2*pi);
    intercept = P(2);
    disp(['intercept is:',num2str(intercept* 180/pi)])
    disp(['fdelta:',num2str(fdelta)])
    disp('Phase offset after mod:')
    disp(phase)
    
    % another way to get phase offset
    avg_phase = mean(phase); 
    disp(['The signal is Phase Shifted with ', num2str(avg_phase * 180/pi), ' degrees'])

    %rx_vec = rx_vec * exp(-1j*avg_phase); % phase syncronization
    rx_vec = rx_vec * exp(-1j*intercept); % phase syncronization
    disp(['The length of rx_vec is', num2str(length(rx_vec))])
    
    % Frequency synchronization
    for i =1:1:length(rx_vec)
        rx_vec(i) = rx_vec(i) * exp(-1j*2*pi*fdelta*Tsamp);
    end
    
    MF_output_cut_without_premable = MF_output_cut_without_premable * exp(-1j*avg_phase); % Match_filter with time sync and phase sync
    rx_vec = rx_vec./abs(rx_vec); % normalise re_vec
    
    disp('Symbols acquired!')
    
    % Minimum Eucledian distance detector
    % Relate the detection to Detection region
    metric = abs(repmat(rx_vec.',1,4) - repmat(const, length(rx_vec), 1)).^2; % compute the distance to each possible symbol
    [tmp, m_hat] = min(metric, [], 2); % find the closest for each received symbol
    m_hat = m_hat'-1;   % get the index of the symbol in the constellation
    m_hat = de2bi(m_hat, 2, 'left-msb')'; %make symbols into bits
    a_hat = m_hat(:)'; %write as a vector
    
    disp('Got the message !!!')
    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
    
    
    bits = a_hat;
    x = rx_vec;
    pulse_train = MF_output_cut_without_premable;
  
    disp(time_value)
    % change the time_value back to 0.2s
    disp('change the time_value back for new preamble detection!')
    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
    % resume(recObj);
    time_value = 0.2;
    set(recObj,'TimerPeriod',time_value);
    disp(['After we got the message we change the time_value back,and the time value now is:',num2str(time_value)])
    recObj.UserData.flag = 0;
    recObj.UserData.time_value = 0.2;
end
          


%------------------------------------------------------------------------------
% HOW TO SAVE DATA FOR THE GUI
%   NOTE THAT THE EXAMPLE HERE IS ONLY USED TO SHOW HOW TO OUTPUT DATA
%------------------------------------------------------------------------------

% Step 1: save the estimated bits
recObj.UserData.pack = bits;

% Step 2: save the sampled symbols
recObj.UserData.const = x;

% Step 3: provide the matched filter output for the eye diagram
recObj.UserData.eyed.r = pulse_train;
recObj.UserData.eyed.fsfd = fsfd;

% Step 4: Compute the PSD and save it. 
% !!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
[pxx, f] = pwelch(pulse_train,1024,768,1024, fs); % note that pwr_spect.f will be normalized frequencies
f = fftshift(f); %shift to be centered around fs
f(1:length(f)/2) = f(1:length(f)/2) - fs; % center to be around zero
p = fftshift(10*log10(pxx/max(pxx))); % shift, normalize and convert PSD to dB
recObj.UserData.pwr_spect.f = f;
recObj.UserData.pwr_spect.p = p;

% In order to make the GUI look at the data, we need to set the
% receive_complete flag equal to 1:
recObj.UserData.receive_complete = 1; 
    
end
