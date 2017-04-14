clear all; close all; clc;

addpath(genpath('../../../../../../tools/mfcc_code/'));
addpath('ola/');

wavpath = '../../wav_48k/';
files = dir(strcat(wavpath,'*.wav'));

tgtdir1 = '../../feats/lmsp/';
tgtdir2 = '../../feats/melsp/';
mkdir(tgtdir1);
mkdir(tgtdir2);

fs = 24000;
n1ms = fs/1000;
frSize = round(50*n1ms);
frShift = round(1*n1ms);
frOvlap = frSize - frShift;
nfft = 2048;
nfftby2 = round(nfft/2 + 1);
nbands = 80;

win = hann(frSize);
minpha_flag = 0;
num_iters = 50;

for i = 1%:length(files)
    
    % [fname,tok] = strtok(files(i).name,'.');
    % fprintf('Processing file %s ... \n',fname);
    
    [y,fss] = audioread(strcat('TheFirebird_32_007','.wav'));
    y = resample(y,fs,fss);
    soundsc(y,fs);
    pause(length(y)/fs)
    %     y = filter([1 -0.97], 1, y);
    %     soundsc(y,fs);
    %     pause(length(y)/fs)
    
    yb = buffer(y,frSize,frOvlap,'nodelay');
    ybw = bsxfun(@times,yb,hann(frSize));
    Y = (fft(ybw,nfft));
    magy = abs(Y);
    phay = unwrap(angle(Y));
    magy = magy + 1e-6;
    hpY = magy(1:nfftby2,:).^2;
    lhY = 20*log10(magy(1:nfftby2,:));
    magyr = 10.^(lhY/20);
    magyr = [magyr;magyr(end-1:-1:2,:)];
    
    yr = real(ifft(magyr.*exp(1i*phay),nfft));
    yr = ola(yr,frShift,@hann,[]);
    soundsc(yr,fs);
    pause(length(y)/fs)
    
    
    [y_ola,y_olaar,y_olagl] = griffin_lim((magyr),frSize,frShift,frOvlap,win,minpha_flag,num_iters);
    
    
    
    soundsc(y_ola/max(abs(y_ola)),fs);
    pause(length(y)/fs)
    soundsc(y_olaar/max(abs(y_olaar)),fs);
    pause(length(y)/fs)
    soundsc(y_olagl/max(abs(y_olagl)),fs);
    
    
    [y_ola,y_olaar,y_olagl] = griffin_lim((magyr.^(1.2)),frSize,frShift,frOvlap,win,minpha_flag,num_iters);
    
    
    
    soundsc(y_ola/max(abs(y_ola)),fs);
    pause(length(y)/fs)
    soundsc(y_olaar/max(abs(y_olaar)),fs);
    pause(length(y)/fs)
    soundsc(y_olagl/max(abs(y_olagl)),fs);
    
    
    
    [y_ola,y_olaar,y_olagl] = griffin_lim((magyr.^(1.4)),frSize,frShift,frOvlap,win,minpha_flag,num_iters);
    
    
    
    soundsc(y_ola/max(abs(y_ola)),fs);
    pause(length(y)/fs)
    soundsc(y_olaar/max(abs(y_olaar)),fs);
    pause(length(y)/fs)
    soundsc(y_olagl/max(abs(y_olagl)),fs);
    
    k = 1; ax(k) = subplot(511);plot(y); axis tight;
    k = k + 1; ax(k) = subplot(512);plot(yr); axis tight;
    k = k + 1; ax(k) = subplot(513);plot(y_ola); axis tight;
    k = k + 1; ax(k) = subplot(514);plot(y_olaar); axis tight;
    k = k + 1; ax(k) = subplot(515);plot(y_olagl); axis tight;
    linkaxes(ax,'x');
    
    %     % mel-spec
    %     melsp = audspec(hpY, fs, nbands, 'mel');
    %     melsp = 10*log10(melsp);
    
    %     dlmwrite(strcat(tgtdir1,fname,'.lmsp'),lhY','delimiter',' ');
    %     dlmwrite(strcat(tgtdir2,fname,'.melsp'),melsp','delimiter',' ');
    
end


