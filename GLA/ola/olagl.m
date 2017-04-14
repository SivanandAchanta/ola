function Y = olagl(X,sShift,hWindow,Dim)
% ola       - Overlap and Add using Griffin & Lim's method
%
%   Y = olagl(X,sShift,hWindow,Dim)
%
%   Performs the overlap and add method on the input matrix X to create the
%   output array Y using Griffin & Lim's method.
%
%   Example
%       Fs = 10000;                 % Sampling frequency
%       t = (0:2*Fs)/Fs;              % Time vector
%
%       A = [0.25 0.5 1 0.5 0.25];  % Amplitudes
%       F = [100 125 150 175 425];  % Frequencies
%       x = A*sin(2*pi*F.'*t);      % Signal
%
%       Framesize = 256;                        % Frame size
%       pOverlap = 0.4;                         % Precentage overlap
%       sOverlap = fix(Framesize*(1-pOverlap)); % Sample overlap
%       sShift = Framesize - sOverlap;          % Sample shift
%       hWindow = @hamming;                     % Windowing function
%
%       X = buffer(x,Framesize,sOverlap,'nodelay');     % Framed signal
%       nFrames = size(X,2);                            % Number of frames
%       X = X.*repmat(hWindow(Framesize),1,nFrames);    % Apply window
%
%       % Normal overlap and add
%       Y1 = ola(X,sShift,[],1);             % Reconstructed signal
%
%       % Allen & Rabiner's method
%       Y2 = ola(X,sShift,@hamming,1);
%       Y2 = olaar(X,sShift,@hamming,1);
%
%       % Griffin & Lim's method
%       Y3 = olagl(X,sShift,@hamming,1);
%
%       figure
%       hold all
%       plot(x,'bx')
%       plot(Y1,'r')
%       plot(Y2,'g')
%       plot(Y3,'c')
%       legend('Original','ola','olaar','olagl')
%
%   See also
%       ola olaar olagl
%

%% Author Information
% Pierce Brady
% Smart Systems Integration Group - SSIG
% Cork Institute of Technology, Ireland.
%

%% Reference
% SSBoll79 function by Esfandiar Zavarehei
% FRAMES2VEC function by Kamil Wojcicki, UTD, July 2011
%

%% Assign defaults
if nargin<4 || isempty(Dim), Dim = 1; end
if nargin<3 || isempty(hWindow), hWindow = []; end

%%
if Dim == 1
    [Framesize, nFrames] = size(X);
    if nargin<2 || isempty(sShift), sShift = Framesize/2; end
    Y = zeros((nFrames-1)*sShift+Framesize,1);
    if isempty(hWindow)
        % No window, so perform normal ola
        for i = 1:nFrames
            % Loop through each frame
            iStart = (i-1)*sShift + 1;
            iFinish = iStart + Framesize - 1;
            Y(iStart:iFinish) = Y(iStart:iFinish) + X(:,i);
        end
    else
        W = zeros((nFrames-1)*sShift+Framesize,1);
        w = feval(hWindow,Framesize);
        X = diag(w)*X;
        w = w.^2;
        for i = 1:nFrames
            % Loop through each frame
            iStart = (i-1)*sShift + 1;
            iFinish = iStart + Framesize - 1;
            Y(iStart:iFinish) = Y(iStart:iFinish) + X(:,i);
            W(iStart:iFinish) = W(iStart:iFinish) + w;
        end
        Y = Y./W;
    end
elseif Dim == 2
    [nFrames,Framesize] = size(X);
    if nargin<2 || isempty(sShift), sShift = Framesize/2; end
    Y = zeros(1,(nFrames-1)*sShift+Framesize);
    if isempty(hWindow)
        % No window, so perform normal ola
        for i = 1:nFrames
            % Loop through each frame
            iStart = (i-1)*sShift + 1;
            iFinish = iStart + Framesize - 1;
            Y(iStart:iFinish) = Y(iStart:iFinish) + X(i,:);
        end
    else
        W = zeros(1,(nFrames-1)*sShift+Framesize);
        w = feval(hWindow,Framesize).';
        X = X*diag(w);
        w = w.^2;
        for i = 1:nFrames
            % Loop through each frame
            iStart = (i-1)*sShift + 1;
            iFinish = iStart + Framesize - 1;
            Y(iStart:iFinish) = Y(iStart:iFinish) + X(i,:);
            W(iStart:iFinish) = W(iStart:iFinish) + w;
        end
        Y = Y./W;
    end
end

end
