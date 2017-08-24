function [rtfr, tfrftic, tfrttic] = sqSTFTCT_fast(x, Bin, Hop, fs, h0, Dh0, MT, varargin)
% Synchrosqueezing by Li Su, 2015
%
% [tfr, rtfr, tfrftic, tfrttic] = sqSTFT_fast(x, 5, Hop, 44100, h, Dh, 30, 'HighFreq', 0.2, 'Smooth', 2, 'Reject', 3)
%
% Input:
% x: input signal (a vector)
% Bin: frequency resolution in Hz (the frequency difference between two consecutive bins)
% Hop: time resolution in second
% fs: sampling frequency
% h: window function
% Dh: 1st order derivative of the window function
% HighFreq: the highest frequency for the output, between 0 and 0.5 (Nyquist rate)
% Smooth: smoothed synchrosqueezed transform (default:0)
% Reject: windowed rejection parameter, 0 = no rejection (default:0)
% Output:
% tfr: STFT (complex valued)
% rtfr: SST (complex valued)
% tfrftic: index of frequency axis
% tfrttic: index of time axis

HighFreq = 0.5;
Smooth = 0;
Reject = 0;
for var_i = 1:length(varargin)
    if strcmp(varargin{var_i}, 'HighFreq')
        HighFreq = varargin{var_i + 1};
    end
    if strcmp(varargin{var_i}, 'Smooth')
        Smooth = varargin{var_i + 1};
    end
    if strcmp(varargin{var_i}, 'Reject')
        Reject = varargin{var_i + 1};
    end
end

if HighFreq > 0.5
    error('Out of Nyquist rate');
end

if Reject < 0
    error('Invalid rejection window width');
end

if Smooth < 0 
    error('Invalid smoothing number');
end

if size(h0,2) > size(h0,1)
    h0 = h0';
end

if size(Dh0,2) > size(Dh0,1)
    Dh0 = Dh0';
end

if size(x,2) > size(x,1)
    x = x';
end

dim = size(h0,2);
	% for tfr
    
alpha = Bin/fs;
Hop = round(fs*Hop);
N = length(-0.5+alpha:alpha:0.5);
Win_length = max(size(h0));
% TH = 5*N/size(h0,1); 
tfrtic = linspace(0, 0.5, round(N/2))' ;
	% for tfrsq
% Lidx = ceil( (N/2)*(lowFreq/0.5) ) ; 
% Hidx = floor( (N/2)*(highFreq/0.5) ) ; 
% fLen = Hidx - Lidx + 1 ;

Overlap = Win_length-Hop;
x_Frame = buffer(x, Win_length, Overlap);
x_Frame = x_Frame(:,ceil(Overlap/Hop/2)+1:end);
Stime = round(Hop-Win_length/2+1)+ceil(Overlap/Hop/2)*Hop;
tfrttic = Stime:Hop:Stime+(size(x_Frame,2)-1)*Hop;

Hidx = floor( (N/2)*(HighFreq/0.5) )+1 ;
tfrftic = linspace(0, HighFreq, Hidx)' ;
K = Hidx; L = size(x_Frame,2);
rtfr = zeros(K, L) ;

% s = rng('shuffle');
for ii = 1:MT
fprintf(['ConceFT total: ',num2str(MT),'; now: %4d\n'], ii) ;

rv = randn(1, dim) ; rv = rv ./ norm(rv) ;
rh = rv * h0';
rDh = rv * Dh0';
h = rh'; Dh = rDh';

tf0 = x_Frame.*repmat(h, [1 L]); 	% for h
tf0 = fft(tf0, N, 1);
tf0 = tf0(1:K, :);

tf2 = x_Frame.*repmat(Dh, [1 L]); 	% for h
tf2 = fft(tf2, N, 1);
tf2 = tf2(1:K, :);

% avoid_warn = find(tf0~=0);
omega = imag(N*tf2./tf0/(2.0*pi));
omega(isinf(omega)|isnan(omega))=0;

if Reject > 0
    TH = Reject*N/size(h,1); 
    omega(abs(omega)>TH/2)=0;
end

omega = round(omega);
% omega(:,:,ii) = reshape(omega_temp,size(tf2,1), size(tf2,2));
% fprintf('\n') ;

% omega(abs(omega)>TH/2)=0;
OrigIndex = repmat((1:K)', [1 L]);
omega(OrigIndex - omega < 1+2*Smooth | OrigIndex - omega > K-2*Smooth)=0;
% ReasIndex = OrigIndex - omega;

Ex = mean(abs(x).^2);
Threshold = 1.0e-8*Ex;	% originally it was 1e-6*Ex
tf0(abs(tf0) < Threshold) = 0;

totLength = size(tf0,1)*size(tf0,2);
new_idx = (1:totLength)'-omega(:);
if Smooth == 0
    rtf0 = accumarray([1; new_idx(2:totLength-1); totLength],tf0(:));
else
    SmoothWin = triang(1+2*Smooth)./sum(triang(1+2*Smooth));
    rtf0 = accumarray([1; new_idx(2:totLength-1); totLength], SmoothWin(1).*tf0(:));
    for jj = 1:Smooth-1
        rtf0 = rtf0 + accumarray([1:jj; new_idx(jj+1:totLength-jj)-jj; totLength-jj+1:totLength], SmoothWin(jj+1).*tf0(:));
        rtf0 = rtf0 + accumarray([1:jj; new_idx(jj+1:totLength-jj)+jj; totLength-jj+1:totLength], SmoothWin(jj+1).*tf0(:));
    end
end
% rtf0 = accumarray((1:totLength)'-omega(:), abs(tf0(:)));
rtf0 = [rtf0; zeros(totLength-length(rtf0),1)];
rtfr = rtfr + reshape(rtf0, K, L);
% rtfr(:,:,ii) = reshape(rtf0, K, L);
end
rtfr = rtfr./MT;

% omega(isinf(omega)|isnan(omega))=0;

% if strcmp(optM,'mean')
%     rtfr = mean(rtfr,3);
% elseif strcmp(optM,'median')
%     rtfr = median(rtfr,3);
% end

end