% Purpose : An example script to demostrate time-scale modification (TSM) using OLA and SOLA

% References
% [1] S. Roucos, A. Wilgus, ‘High quality Time-Scale Modification of Speech’, ICASSP-85, pp. 236-239, 1985.

clear all; close all; clc;

% Analysis conditions
fs        = 16000;
frSizems  = 30;
frShiftms = 10;
frSize    = (frSizems/1000)*fs;
frShift   = (frShiftms/1000)*fs;
frOvlap   = frSize - frShift;
nfft      = 512;
nfftby2   = round(nfft/2 + 1);
hfhz      = linspace(0,fs/2,nfftby2);
w         = hamming(frSize);
max_delay = 80; % SOLA parmaeter

% Plot flags
pflag = 1;
plot_crosscorr = 1;

% Matlab plot parameters
font_size = 14;

% Time scale modification factor
tsmf = 0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a synthetic signal (pulse train or a rectangular wave)
t_on      = 10;
t_off     = 160 - t_on;
y = [zeros(1,40) ones(1,t_on) zeros(1,t_off) ones(1,t_on) zeros(1,t_off) ones(1,t_on) zeros(1,t_off) ones(1,t_on) zeros(1,t_off) ones(1,t_on) zeros(1,t_off) ones(1,t_on) zeros(1,t_off) ones(1,t_on) zeros(1,t_off) ones(1,t_on) zeros(1,t_off)];
y = y(:);
y = [y;y(41:end)];


% Call OLA
[ysn_ola] = ola_v3(y,frSize,frShift,w,tsmf);

% Call SOLA
[ysn_sola] = sola(y,frSize,frShift,w,tsmf,max_delay,plot_crosscorr);

% Plot the result
if pflag
    ax(1) = subplot(311); plot((1:length(y))/fs,y,'linewidth',1.2); axis tight; title('Original signal','FontSize',font_size,'FontWeight','bold');
    ax(2) = subplot(312); plot((1:length(ysn_ola))/fs,ysn_ola,'r-','linewidth',1.2); axis tight; title(['Time scale modified signal using OLA method with stretch factor = ', num2str(tsmf)],'FontSize',font_size,'FontWeight','bold');
    ax(3) = subplot(313); plot((1:length(ysn_sola))/fs,ysn_sola,'m-','linewidth',1.2); axis tight; title(['Time scale modified signal using SOLA method with stretch factor = ', num2str(tsmf)],'FontSize',font_size,'FontWeight','bold');
    xlabel('time(sec)','FontSize',font_size,'FontWeight','bold');
    set(ax,'fontweight','bold');
end
