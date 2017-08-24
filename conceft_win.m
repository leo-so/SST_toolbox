% Generating an overcomplete set of windows from the Gaussian-Hermite series
% N: total number of tapers
% order: number of the orthogonal window basis
% H: the set of windows for ConceFT

function H = conceft_win(N, order, MT)
% h=h./norm(h);
% for oi = 2:order
%     dh = dwindow(h(:,end));
%     dh = dh./norm(dh);
%     h = [h dh];
% end
[h, ~, ~] = hermf(N, order, 3);
h=h';

H = [];

for n = 1:MT
    rv = randn(1, order) ; 
    rv = rv ./ norm(rv) 
    rh = rv * h';
    H = [H rh'];
end