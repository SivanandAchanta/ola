function [y,y1,y2] = griffin_lim(X,frSize,frShift,frOvlap,win,minpha_flag,num_iters)


[nfft,nof] = size(X);
len_sig = frShift*nof;


if minpha_flag
    nfftby2 = round(nfft/2 + 1);
    C = real(ifft(log(X),nfft));
    
    % assymetrize cestrum
    C(1:nfftby2,:) = 2*C(1:nfftby2,:);
    C(nfftby2+1:end,:) = 0;
    
    iC = real(ifft(exp(fft(C,nfft)),nfft));
    iC = iC(1:frSize,:);
    y = ola(iC,frShift,@hann,[]);
    y = y(1:len_sig);
    y1 = y;
    y2 = y;
else
    y = randn(len_sig,1);
    y1 = y;
    y2 = y;
end


for i = 1:num_iters
    
    yw = buffer(y,frSize,frOvlap);    
    yw = bsxfun(@times,yw,win);
    Y = fft(yw,nfft);
    phay = unwrap(angle(Y));    
    MY = X.*exp(1i*phay);
    my = real(ifft(MY,nfft));
    my = my(1:frSize,:);    
    
    yw = buffer(y1,frSize,frOvlap);    
    yw = bsxfun(@times,yw,win);
    Y = fft(yw,nfft);
    phay = unwrap(angle(Y));    
    MY = X.*exp(1i*phay);
    my1 = real(ifft(MY,nfft));
    my1 = my1(1:frSize,:);    
    
    yw = buffer(y2,frSize,frOvlap);    
    yw = bsxfun(@times,yw,win);
    Y = fft(yw,nfft);
    phay = unwrap(angle(Y));    
    MY = X.*exp(1i*phay);
    my2 = real(ifft(MY,nfft));
    my2 = my2(1:frSize,:);    
        
    y = ola(my,frShift,@hann,[]); 
    y1 = olaar(my1,frShift,@hann,[]); 
    y2 = olagl(my2,frShift,@hann,[]); 
    
    y = y(240:len_sig+240-1);
    y1 = y1(240:len_sig+240-1);
    y2 = y2(240:len_sig+240-1);
    

    
end


