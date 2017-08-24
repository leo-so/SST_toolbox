[x, fs] = wavread('C:\Users\user\Dropbox\中研院合作音樂家003\DATA\Symphony\SY05_Berlioz\SY05_Berlioz_audio.wav');
x = mean(x,2);
y = x(1*fs:fs*11);
% [h,Dh,tt] = hermf(4097, 3, 6);
h = hann(4097)';
Dh = dwindow(h')';
[tfr, rtfr, tfrftic, tfrttic] = sqSTFT_fast(y, 5, 0.01, fs, h(1,:), Dh(1,:), 'HighFreq', 0.2, 'Smooth', 0, 'Reject', 0);
% [tfr, rtfr2, tfrftic, tfrttic] = rsSTFT_fast(x, 1, 0.01, fs, h(1,:), Dh(1,:), 'HighFreq', 0.1, 'Smooth', 0, 'Reject', 0);
% [conceft, tfrftic, tfrttic] = sqSTFTCT_fast(x, 1, 0.01, fs, h, Dh, 30, 'HighFreq', 0.1, 'Smooth', 0, 'Reject', 0);
% rtfr(rtfr==0)=1E-3;
figure(1)
imagescQ(tfr, tfrftic, tfrttic, fs); axis xy;
% title('STFT');
xlabel('time(s)'); ylabel('frequency (Hz)');
% wavwrite(y,fs,'sym.wav');
% figure(2)
% imagescQ(rtfr, tfrftic, tfrttic, fs);
% title('SST');
% figure(3)
% imagescQ(rtfr2, tfrftic, tfrttic, fs);
% title('RM');