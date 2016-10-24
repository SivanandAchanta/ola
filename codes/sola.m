function [ysn_sola] = sola(y,frSize,frShift,win,tsmf,max_delay,plot_crosscorr)

% Purpose : To implement SOLA method for time-scale modification fo speech
% signals

% Inputs:
% [1] y - original signal
% [2] frSize - frame size in samples (usually 30ms*sr)
% [3] frShift - frame shift in samples (usually 10ms*sr)
% [4] win - hamming window of length frame size
% [5] tsmf - time-scale modification factor (0.5 to 1.3)
% [6] max_delay - maximum allowable delay for cross-correlation
% [7] plot_crosscorr - Plotting and visualizing the details of SOLA steps 

% Outputs:
% [1] ysn_sola - sola modified signal

% References
% [1] S. Roucos, A. Wilgus, ‘High quality Time-Scale Modification of Speech’, ICASSP-85, pp. 236-239, 1985.


% Author : Sivanand Achanta


% Set analysis and synthesis frame anchor points
len_synsig = ceil(tsmf*length(y));
ati = 0:frShift:length(y)-frSize; % analysis time instants
sti = 0:(tsmf*frShift):(len_synsig)-frSize; % synthesis time instants
v = win.^2;

% Initialization of Sola signal and normalization signal
ys_sola = zeros(len_synsig,1);
normalization_sig = zeros(len_synsig,1) + 1e-6;
delay = 0;

if plot_crosscorr
    h1 = figure(1);
end

for i = 1:length(sti)
    
    sp_nodelay = sti(i) + 1;
    ep_nodelay = sti(i) + frSize;
    
    if (i > 1)
        
        % Take the windowed analysis segment
        x_w = win.*y(ati(i)+1:ati(i)+frSize);
        
        % Do the cross-correlation in both directions (Actually one of them is redundant will remove in the next revision)
        r1 = xcorr(x_w,ys_sola(sp_nodelay:ep_nodelay));
        r2 = xcorr(ys_sola(sp_nodelay:ep_nodelay),x_w);
        
        % Compute the maximum correlation and delay in each of the direction
        [maxval1,d1] = max(r1(frSize:frSize+max_delay));
        d1 = d1 - 1;
        [maxval2,d2] = max(r2(frSize:frSize+max_delay));
        d2 = d2 - 1;
        
        % Take the delay which gives maximum correlation
        if (maxval1 > maxval2)
            delay = -d1; % note the minus sign here ... ***
        else
            delay = d2;
        end
        
        % Plot
        if plot_crosscorr
            nr = 6; nc = 1;
            
            % Plot the signal along with the current analysis window
            k = 1; ax(k) = subplot(nr,nc,k); plot((1:ati(i)+frSize),y(1:ati(i)+frSize)); hold on; plot((ati(i)+1:ati(i)+frSize),win,'r-'); hold off; axis tight;
            title('Signal along with the current analysis window','FontSize',12,'FontWeight','bold');
            
            % Place the analysis segment at the synthesis time instant
            k = k + 1; ax(k) = subplot(nr,nc,k); plot([zeros(sti(i),1);x_w]); axis tight;
            title('Analysis segment placed at the synthesis time instant','FontSize',12,'FontWeight','bold');
            
            % Compensate the delay required for synchronization
            k = k + 1; ax(k) =  subplot(nr,nc,k); plot([zeros(sti(i)+delay,1);x_w]);
            title('Delay compensated synthesis time instant','FontSize',12,'FontWeight','bold');
            
            % Plot the cross correlation fucntion in both directions
            k = k + 1; ax(k) =  subplot(nr,nc,k); plot(r1); axis tight;
            title('Cross-correlation in one direction','FontSize',12,'FontWeight','bold');
            
            k = k + 1; ax(k) =  subplot(nr,nc,k); plot(r2); axis tight;
            title('Cross-correlation in other direction','FontSize',12,'FontWeight','bold');
            
            % Plot the SOLA signal constructed so far
            k = k + 1; ax(k) =  subplot(nr,nc,k); plot(ys_sola(1:ep_nodelay)); axis tight;
            title('SOLA signal constructed so far','FontSize',12,'FontWeight','bold');
            
            linkaxes(ax(1:3),'x');
            xlabel('Samples','FontSize',12,'FontWeight','bold');
            pause
        end
        
    else
        x_w = win.*y(ati(i)+1:ati(i)+frSize);
    end
    
    % Add the delay (here SOLA differs from OLA method)
    sp = sti(i) + 1 + delay ;
    ep = sti(i) + frSize + delay ;
    
    % OLA after synchronization (same as OLA)
    ys_sola(sp:ep) = ys_sola(sp:ep) + win.*x_w;
    normalization_sig(sp:ep) = normalization_sig(sp:ep) + v;
end

ysn_sola = ys_sola./normalization_sig;
