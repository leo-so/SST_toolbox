N = 4097; % window size, odd number is suggsted 
order = 3; % number of Gaussian-Hermite series
MT = 10; % total number of tapers (windows)
H = conceft_win(N, order, MT);

% Then compute the STFTs using H(:,1) to H(:,end)