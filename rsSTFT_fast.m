function [tfr, rtfr, tfrftic, tfrttic] = rsSTFT_fast(x, Bin, Hop, fs, h, Dh, varargin)
% Synchrosqueezed transform (SST) by Li Su, 2016/01/11
%
% [tfr, rtfr, tfrftic, tfrttic] = sqSTFT_fast(x, 5, Hop, 44100, h, Dh, 'HighFreq', 0.2, 'Smooth', 2, 'Reject', 3)
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

if size(h,2) > size(h,1)
    h = h';
end

if size(Dh,2) > size(Dh,1)
    Dh = Dh';
end

if size(x,2) > size(x,1)
    x = x';
end

	% for tfr
alpha = Bin/fs;
Hop = round(fs*Hop);
N = length(-0.5+alpha:alpha:0.5);
Win_length = max(size(h));

	% for tfrsq
% Lidx = ceil( (N/2)*(lowFreq/0.5) ) ; 
% Hidx = floor( (N/2)*(highFreq/0.5) ) ; 
% fLen = Hidx - Lidx + 1 ;

Overlap = Win_length-Hop;
x_Frame = buffer(x, Win_length, Overlap);
x_Frame = x_Frame(:,ceil(Overlap/Hop/2)+1:end);
Stime = round(Hop-Win_length/2+1)+ceil(Overlap/Hop/2)*Hop;
tfrttic = Stime:Hop:Stime+(size(x_Frame,2)-1)*Hop;

tfr = x_Frame.*repmat(h, [1 size(x_Frame,2)]); 	% for h
tfr = fft(tfr, N, 1);
% tfr = tfr(1:round(N/2),:);

tf2 = x_Frame.*repmat(Dh, [1 size(x_Frame,2)]); 	% for h
tf2 = fft(tf2, N, 1);
% tf2 = tf2(1:round(N/2),:);

Th=h.*(-floor(length(h)/2):floor(length(h)/2))';  % for Th
tf3 = x_Frame.*repmat(Th, [1 size(x_Frame,2)]); 	% for h
tf3 = fft(tf3, N, 1);

Hidx = floor( (N/2)*(HighFreq/0.5) )+1 ; 
tfr = tfr(1:Hidx,:); tfr0 = tfr;
tf2 = tf2(1:Hidx,:);
tf3 = tf3(1:Hidx,:);
tfrftic = linspace(0, HighFreq, Hidx)' ;

		% get the first order omega
omega = zeros(size(tf2)) ;
gdf = zeros(size(tf2)) ;
avoid_warn = find(tfr~=0);
omega(avoid_warn) = imag(N*tf2(avoid_warn)./tfr(avoid_warn)/(2.0*pi));
gdf(avoid_warn) = real(tf3(avoid_warn)./tfr(avoid_warn)/Hop);

if Reject > 0
    TH = Reject*N/size(h,1); 
    tfr(abs(omega)>TH/2)=0;
end
tfr(abs(gdf)>(length(h)/Hop/2))=0;

omega = round(omega);
gdf = round(gdf);
OrigIndex = repmat((1:Hidx)', [1 size(tfr,2)]);
omega(OrigIndex - omega < 1+2*Smooth | OrigIndex - omega > Hidx-2*Smooth)=0;
OrigIndex = repmat((1:size(tfr,2)), [size(tfr,1) 1]);
gdf(OrigIndex + gdf < 1 | OrigIndex + gdf > size(tfr,2))=0;
% ReasIndex = OrigIndex - omega;

Ex = mean(abs(x).^2);
Threshold = 1.0e-6*Ex;	% originally it was 1e-6*Ex
tfr(abs(tfr) < Threshold) = 0;

totLength = size(tfr,1)*size(tfr,2);
new_idx = (1:totLength)' + gdf(:)*Hidx -omega(:);
if Smooth == 0
    rtfr = accumarray([1; new_idx(2:totLength-1); totLength],tfr(:));
else
    SmoothWin = triang(1+2*Smooth)./sum(triang(1+2*Smooth));
    rtfr = accumarray([1; new_idx(2:totLength-1); totLength], SmoothWin(1).*tfr(:));
    for ii = 1:Smooth-1
        rtfr = rtfr + accumarray([1:ii; new_idx(ii+1:totLength-ii)-ii; totLength-ii+1:totLength], SmoothWin(ii+1).*tfr(:));
        rtfr = rtfr + accumarray([1:ii; new_idx(ii+1:totLength-ii)+ii; totLength-ii+1:totLength], SmoothWin(ii+1).*tfr(:));
    end
end
rtfr = [rtfr; zeros(totLength-length(rtfr),1)];

rtfr = reshape(rtfr,size(tfr,1),size(tfr,2));
tfr = tfr0;
end