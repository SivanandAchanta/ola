function [ysn_ola] = ola_v3(y,frSize,frShift,win,tsmf)

% Purpose : Perform OLA in time domain

% Inputs: 
% [1] y - original signal
% [2] frSize - frame size in samples (usually 30ms*sr)
% [3] frShift - frame shift in samples (usually 10ms*sr)
% [4] win - hamming window of length frame size
% [5] tsmf - time-scale modification factor

% Outputs:
% [1] ysn_ola - ola modified signal

% References
% [1]: Signal Estimation from Modified Short-Time Fourier Transform --
%       DANIEL W. GRIFFIN AND JAE S. LIM (IEEE ASSP - 1984)

% Author : Sivanand Achanta

% Set analysis and synthesis frame anchor points
len_synsig = ceil(tsmf*length(y));
ati = 0:frShift:length(y)-frSize; % analysis time instants
sti = 0:(tsmf*frShift):(len_synsig)-frSize; % synthesis time instants
v = win.^2;

% Initializing the synthesis signal and the normalization signal
ys_ola = zeros(len_synsig,1);
normalization_sig = zeros(len_synsig,1) + 1e-6;

for i = 1:length(sti)
    % begin and end points of synthesis signal
    sp = sti(i) + 1;
    ep = sti(i) + frSize;
    
    % overlap-add
    ys_ola(sp:ep) = ys_ola(sp:ep) + v.*y(ati(i)+1:ati(i)+frSize);
    
    % normalization signal
    normalization_sig(sp:ep) = normalization_sig(sp:ep) + v;
end
ysn_ola = ys_ola./normalization_sig;
