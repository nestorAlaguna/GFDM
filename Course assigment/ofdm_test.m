close all; clear all; clc;
% SUBCARRIERS:
K = 512; % Number of samples per sub-symbol
Kon = 151; % Number of allocated subcarriers
Kset = 1: Kon; % Vector of allocated subcarriers

% SUB-SYMBOLS:
M = 15; % Number of sub-symbols available
Mon = 14; % Number of sub-symbols allocated
Mset= 2: Mon+1; % Vector of sub-symbols allocated

% Number of blocks, Block length, CP, pulse shape, roll factor
blocks = 10; % Number of blocks used by the signal
blockLength = M*K; % Total allocated space for the signal
CP = ceil(0.1*blockLength);
% CP=0;
a = 0.1; % rolloff factor of the pulse shaping filter
mu = 4; % Number of bits in every QAM symbol
sOFDM = zeros((blocks * blockLength)+CP, 1); % GFDM signal
symbols_tx = [];
symbols_qam_tx = [];

% ser = zeros(1,21); % SER array
biterr_Nr = zeros(1,21); % Number of errored bits
biterr_rate = zeros(1,21); % BER array
snr_max = 20;
snr_interval = 0:snr_max;


[ofdm_s_tx,ofdm_qam_tx,ofdm_s_rx,ofdm_qam_rx,sOFDM,serofdm]= modemod(2,1,1,0.01,mu,1,0);

function [symbols_tx,symbols_qam_tx,symbols_rx,symbols_qam_rx,sGFDM,ser] = modemod(Kop,M,Mon,a,mu,Mset,CP)
%% Creating GFDM waveform paramters

% SUBCARRIERS:
K = 512; % Number of samples per sub-symbol
Kon = 151; % Number of allocated subcarriers
Kset = 1: Kon; % Vector of allocated subcarriers
% Kop = 2; % Number of overlapping subcarriers
% 
% % SUB-SYMBOLS:
% M = 15; % Number of sub-symbols available
% Mon = 14; % Number of sub-symbols allocated
% Mset= 2: Mon+1; % Vector of sub-symbols allocated

% Number of blocks, Block length, CP, pulse shape, roll factor
blocks = 10; % Number of blocks used by the signal
blockLength = M*K; % Total allocated space for the signal
CP = ceil(0.1*blockLength);
% CP=0;
% a = 0.1; % rolloff factor of the pulse shaping filter
mu = 4; % Number of bits in every QAM symbol


%% Allocation of arrays
sGFDM = zeros((blocks * blockLength)+CP, 1); % GFDM signal
symbols_tx = [];
symbols_qam_tx = [];

% ser = zeros(1,21); % SER array
biterr_Nr = zeros(1,21); % Number of errored bits
biterr_rate = zeros(1,21); % BER array
snr_max = 20;
snr_interval = 0:snr_max;

%% MODULATION of GFDM signal with random symbols
%  For every value of SNR in range we calculate and plot the results
for snr_i = 0: snr_max 
    symbols_tx = [];
    symbols_qam_tx = [];
    ser_i = 0;
    ber_i = 0;
    
    % Iterating on the number of blocks
    for block_i = 1:blocks  
        %% Generating random symbols in the interval 0...2 ^ gfdm.mu -1
        rand_symbols = randi(2 ^ mu, Mon * Kon, 1) - 1;
        
        % Applying QAM modulation on the symbols and normalization
        qm2 = qammod(rand_symbols, 2^mu);
        qm2 = qm2 / sqrt(2/3 * (2^mu)); % Normalization
        
        % Adding on every iteration of block the new symbols
        symbols_tx = [symbols_tx; rand_symbols];
        symbols_qam_tx = [symbols_qam_tx; qm2];
        
        % Mapping the linear symbol stream to the GFDM Data matrix and
        % reshaping it according to allocated dimensions
        % Upsampling
        if Kon== K && Mon == M
            D = reshape(qm2, K, M);
        else
            Dm = reshape(qm2, Kon, Mon);
            res1 = zeros(K, Mon);
            res1(Kset, :) = Dm;
            res = zeros(K, M);
            res(:, Mset) = res1;
            D = res;
        end
        
        %% RC - Raised Cosine filter (Time domain) at TX
        %  Defining the filter
        if M == 1
        g = rc_td(M,K);   
        else
        g = rc(M,K,a);
        end
    
        % Applying the modulation in time and convolution
        DD = repmat(K*ifft(D), M, 1);
        x_value = zeros(K*M, 1);
        for m=1:M
            symbol = DD(:, m) .* g;
            symbol = circshift(symbol, K*(m-1));
            x_value = x_value + symbol;
        end
        xcp = [x_value(end - CP + (1:CP),:);x_value];
        %%%%%%%%%%%%% add cpppp
        % Adding the final values to the signal
        sGFDM((block_i-1)*(blockLength+CP)+(1:blockLength+CP)) = xcp;
    end
    
    %% Adding Channel Noise and Demodulation of the signal
    AWGN_sGFDM = awgn(1.*sGFDM, snr_i);
    
    % Demodulation: init of arrays of symbols received
    symbols_rx = [];
    symbols_qam_rx = [];
    
    % Loop on every block
    for block_i = 1:blocks 
       %removing cp
       
        x_block = AWGN_sGFDM((block_i-1)*(blockLength+CP) + (1:(blockLength+CP)));
        rcp = x_block(CP + 1:end);
        % RC - Raised Cosine filter (Frequency domain) at RX
        if M == 1
        g = rc_td(M,K);   
        else
        g = rc(M,K,a);
        end

        g = g(round(1:K/Kop:end));
        G = conj(fft(g));
   
        L = length(G) / M;
        Xhat = fft(rcp);
        Dhat = zeros(K, M);
        for k=1:K
            sc = circshift(Xhat, ceil(L*M/2) - M*(k-1));
            sc = fftshift(sc(1:L*M));
            sc_Matched = sc .* G;
            dhat = ifft(sum(reshape(sc_Matched, M, L), 2)/L);
            Dhat(k,:) = dhat;
        end

        % Downsampling
        if Kon == K && Mon == M
            recv_symbols = reshape(Dhat, numel(Dhat), 1);
        else
            Dm = Dhat(Kset, Mset);
            recv_symbols = reshape(Dm, numel(Dm), 1);
        end
    
    % Unmapping the Dhat matrix to symbols stream
    recv_symbols = recv_symbols * sqrt(2/3 * (2^mu - 1));
    shm = qamdemod(recv_symbols, 2^mu);
    
    % Adding on every loop of block the new symbols
    symbols_rx = [symbols_rx; shm];
    symbols_qam_rx = [symbols_qam_rx; recv_symbols];
    end
    
    % Loop on every symbol to check for symbol errors
    for i = 1 : length(symbols_rx)
        if symbols_rx(i) ~= symbols_tx(i)
            ser_i = ser_i + 1;      
        end
    end
    
    ser(1,snr_i+1) = ser_i/length(symbols_rx);
    [biterr_Nr(1,snr_i+1), biterr_rate(1,snr_i+1)] = biterr(symbols_tx,symbols_rx);

    
end

end

function g = rc(M, K, a)
% RC - Return Raised Cosine filter (time domain)
%
t = linspace(-M/2, M/2, M*K+1); t = t(1:end-1); t = t';

g = (sinc(t) .* cos(pi*a*t) ./ (1-4*a*a*t.*t));
g = fftshift(g);
g(K+1:K:end) = 0;

g = g / sqrt(sum(g.*g));
end

function g = rc_td(M,K)
      
g = zeros(M*K, 1);
g(1:K) = 1;
g = g / sqrt(sum(g.*g));
        
end