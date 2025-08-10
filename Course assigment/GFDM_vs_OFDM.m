%% Creating GFDM waveform paramters
%  And storing all of them into a struct variable

gfdm = struct;
% SUBCARRIERS:
gfdm.K = 512; % Number of samples per sub-symbol
gfdm.Kon = 151; % Number of allocated subcarriers
gfdm.Kset = 1: gfdm.Kon; % Vector of allocated subcarriers
gfdm.Kop = 5; % Number of overlapping subcarriers

% SUB-SYMBOLS:
gfdm.M = 15; % Number of sub-symbols available
gfdm.Mon = 14; % Number of sub-symbols allocated
gfdm.Mset= 1: gfdm.Mon; % Vector of sub-symbols allocated

% Number of blocks, Block length, CP, pulse shape, roll factor
gfdm.blocks = 10; % Number of blocks used by the signal
gfdm.blockLength = gfdm.M*gfdm.K; % Total allocated space for the signal
gfdm.CPLength = ceil(0.01*gfdm.blockLength);
gfdm.a = 0.1; % rolloff factor of the pulse shaping filter
gfdm.mu = 4; % Number of bits in every QAM symbol:
% e.g. gfdm.mu = 4 => 4 bits/symbol => 16-QAM modulation
% e.g. gfdm.mu = 3 => 3 bits/symbol => 8-QAM modulation


%% Creating OFDM waveform paramters
%  And storing all of them into a struct variable

ofdm = struct;
% SUBCARRIERS:
ofdm.K = 512; % Number of subcarriers
ofdm.Kon = 151; % Number of allocated subcarriers
ofdm.Kset = 1:ofdm.Kon; % The vector of allocated subcarriers
ofdm.Kop = 2; % Number of overlapped subcarriers

% In OFDM we do not have sub-symbols per subcarrier
ofdm.M = 1;
ofdm.Mon = ofdm.M;

% Number of blocks, Block length, CP, pulse shape, roll factor
ofdm.blocks = 10; % Number of blocks used by the signal
ofdm.blockLength = gfdm.blockLength; % Total allocated space for the signal
ofdm.CPLength = 0;% ceil(0.1*ofdm.blockLength);
ofdm.pulse = 'rc'; % Raised cosine shaping filter
ofdm.a = 0.1; % rolloff factor of the pulse shaping filter
ofdm.mu = 4; % Number of bits in every QAM symbol:

%% Allocation of arrays
sGFDM = zeros(gfdm.blocks * gfdm.blockLength + gfdm.CPLength, 1); % GFDM signal
sOFDM = zeros(ofdm.blocks * ofdm.blockLength + ofdm.CPLength, 1); % OFDM signal
symbols_gfdm_tx = [];
symbols_qam_tx = [];

ser = zeros(1,21); % SER array
biterr_Nr = zeros(1,21); % Number of erred bits
biterr_rate = zeros(1,21); % BER array
snr_max = 20;
snr_interval = 0:snr_max;

%% MODULATION of GFDM signal with random symbols
%  For every value of SNR in range we calculate and plot the results
for snr_i = 0: snr_max 
    symbols_gfdm_tx = [];
    symbols_qam_tx = [];
    ser_i = 0;
    ber_i = 0;
    
    % Iterating on the number of blocks
    for block_i = 1:gfdm.blocks  
        %% Generating random symbols in the interval 0...2 ^ gfdm.mu -1
        rand_symbols = randi(2 ^ gfdm.mu, gfdm.Mon * gfdm.Kon, 1) - 1;
        
        % Applying QAM modulation on the symbols and normalization
        qm2 = qammod(rand_symbols, 2^gfdm.mu);
        qm2 = qm2 / sqrt(2/3 * (2^gfdm.mu)); % Normalization
        
        % Adding on every iteration of block the new symbols
        symbols_gfdm_tx = [symbols_gfdm_tx; rand_symbols];
        symbols_qam_tx = [symbols_qam_tx; qm2];
        
        % Mapping the linear symbol stream to the GFDM Data matrix and
        % reshaping it according to allocated dimensions
        % Upsampling
        if gfdm.Kon== gfdm.K && gfdm.Mon == gfdm.M
            D = reshape(qm2, gfdm.K, gfdm.M);
        else
            Dm = reshape(qm2, gfdm.Kon, gfdm.Mon);
            res1 = zeros(gfdm.K, gfdm.Mon);
            res1(gfdm.Kset, :) = Dm;
            res = zeros(gfdm.K, gfdm.M);
            res(:, gfdm.Mset) = res1;
            D = res;
        end
        
        %% RC - Raised Cosine filter (Time domain) at TX
        %  Defining the filter
        t = linspace(-gfdm.M/2, gfdm.M/2, gfdm.M*gfdm.K + 1); 
        t = t(1:end-1); 
        t = t';
        g = (sinc(t) .* cos(pi*gfdm.a*t) ./ (1-4*gfdm.a*gfdm.a*t.*t));
        g = fftshift(g);
        g(gfdm.K+1 : gfdm.K : end) = 0;
        g = g / sqrt(sum(g.*g));
    
        % Applying the modulation in time and convolution
        DD = repmat(gfdm.K*ifft(D), gfdm.M, 1);
        x_value = zeros(gfdm.K*gfdm.M, 1);
        for m=1:gfdm.M
            symbol = DD(:, m) .* g;
            symbol = circshift(symbol, gfdm.K*(m-1));
            x_value = x_value + symbol;
        end
    %%%%%%%%%%%%% add cpppp
        % Adding the final values to the signal
        xcp = [x_value(end - gfdm.CPLength + (1:gfdm.CPLength),:);x_value];
        sGFDM((block_i-1)*(gfdm.blockLength + gfdm.CPLength)+(1:(gfdm.blockLength + gfdm.CPLength))) = xcp;
    end
    
    %% Adding Channel Noise and Demodulation of the signal
    AWGN_sGFDM = awgn(1.*sGFDM, snr_i);
    
    % Demodulation: init of arrays of symbols received
    symbols_gfdm_rx = [];
    symbols_qam_rx = [];
    
    % Loop on every block
    for block_i = 1:gfdm.blocks 
        x_block = AWGN_sGFDM((block_i-1)*(gfdm.blockLength + gfdm.CPLength) + (1:(gfdm.blockLength + gfdm.CPLength)));
        rcp = x_block(gfdm.CPLength + 1:end);
        % RC - Raised Cosine filter (Frequency domain) at RX
        t = linspace(-gfdm.M/2, gfdm.M/2, gfdm.M*gfdm.K+1); 
        t = t(1:end-1); 
        t = t';
        g = (sinc(t) .* cos(pi*gfdm.a*t) ./ (1-4*gfdm.a*gfdm.a*t.*t));
        g = fftshift(g);
        g(gfdm.K + 1:gfdm.K:end) = 0;
        g = g / sqrt(sum(g.*g));

        g = g(round(1:gfdm.K/gfdm.Kop:end));
        G = conj(fft(g));
   
        L = length(G) / gfdm.M;
        Xhat = fft(rcp);
        Dhat = zeros(gfdm.K, gfdm.M);
        for k=1:gfdm.K
            sc = circshift(Xhat, ceil(L*gfdm.M/2) - gfdm.M*(k-1));
            sc = fftshift(sc(1:L*gfdm.M));
            sc_Matched = sc .* G;
            dhat = ifft(sum(reshape(sc_Matched, gfdm.M, L), 2)/L);
            Dhat(k,:) = dhat;
        end

        % Downsampling
        if gfdm.Kon == gfdm.K && gfdm.Mon == gfdm.M
            recv_symbols = reshape(Dhat, numel(Dhat), 1);
        else
            Dm = Dhat(gfdm.Kset, gfdm.Mset);
            recv_symbols = reshape(Dm, numel(Dm), 1);
        end
    
    % Unmapping the Dhat matrix to symbols stream
    recv_symbols = recv_symbols * sqrt(2/3 * (2^gfdm.mu - 1));
    shm = qamdemod(recv_symbols, 2^gfdm.mu);
    
    % Adding on every loop of block the new symbols
    symbols_gfdm_rx = [symbols_gfdm_rx; shm];
    symbols_qam_rx = [symbols_qam_rx; recv_symbols];
    end
    
    % Loop on every symbol to check for symbol errors
    for i = 1 : length(symbols_gfdm_rx)
        if symbols_gfdm_rx(i) ~= symbols_gfdm_tx(i)
            ser_i = ser_i + 1;      
        end
    end
    
    ser(1,snr_i+1) = ser_i/length(symbols_gfdm_rx);
    [biterr_Nr(1,snr_i+1), biterr_rate(1,snr_i+1)] = biterr(symbols_gfdm_tx,symbols_gfdm_rx);
    
end

%% Eye diagram plot, scatter plot of the symbols
    eyediagram(symbols_qam_rx(1:100), 2); % figure(1)
    scatterplot(symbols_qam_rx); % figure(2)
    title('Received Signal');
    scatterplot(symbols_qam_tx); % figure(3)
    title('Transmitted Signal');

%% Plotting SER of a GFDM signal
%  Comparing the actual results with the theoretical ones,
%  computed with the formula of erfc
    figure(4);
    ser_theory = 3/2*erfc(sqrt(0.1*(10.^(snr_interval/10))));
    semilogy(snr_interval, ser_theory,'-x');
    hold on;
    semilogy(snr_interval, ser, 'o');
    legend('theory','experiment');
    xlabel('SNR, dB');
    ylabel('Symbol Error Rate');
    hold off;

%% Plotting BER of a GFDM signal 
    figure(5);
    plot(snr_interval, biterr_rate);
    xlabel('SNR, dB');
    title('Bit Error Rate');

%% Create an OFDM signal
% sOFDM = zeros(size(sGFDM));
% for b = 1 : ofdm.blocks 
%     for m = 1 : gfdm.M  
%         rand_symbols_2 = randi(2 ^ ofdm.mu, ofdm.Mon * ofdm.Kon, 1) - 1;
%         qm_2 = qammod(rand_symbols_2, 2^ofdm.mu)/sqrt(2/3*(2^ofdm.mu - 1));
%         
%         % Upsampling
%         if ofdm.Kon == ofdm.K
%             D_2 = reshape(qm_2, ofdm.K, ofdm.M);
%         else
%             Dm_2 = reshape(qm_2, ofdm.Kon, ofdm.Mon);
%             res1 = zeros(ofdm.K, ofdm.Mon);
%             res1(ofdm.Kset, :) = Dm_2;
%             res = zeros(ofdm.K, ofdm.M);
%             res(:, ofdm.Mon) = res1;
%             D_2 = res;
%         end
% 
%         %% RC - Raised Cosine filter (Time domain) at TX
%         %  Defining the filter
%         t = linspace(-ofdm.M/2, ofdm.M/2, ofdm.M*ofdm.K + 1); 
%         t = t(1:end-1); 
%         t = t';
%         g_2 = (sinc(t) .* cos(pi*ofdm.a*t) ./ (1-4*ofdm.a*ofdm.a*t.*t));
%         g_2 = fftshift(g_2);
%         g_2(ofdm.K+1:ofdm.K:end) = 0;
%         g_2 = g_2 / sqrt(sum(g_2.*g_2));
%         
%         % Applying the modulation in time and convolution
%         DD_2 = repmat(ofdm.K*ifft(D_2), ofdm.M, 1);
%         x2 = zeros(ofdm.K*ofdm.M, 1);
%         for m=1:ofdm.M
%             symbol = DD_2(:, m) .* g_2;
%             symbol = circshift(symbol, ofdm.K*(m-1));
%             x2 = x2 + symbol;
%         end
%      
%          xcp1 = [x2(end - ofdm.CPLength + (1:ofdm.CPLength),:);x2];
%         sOFDM((block_i-1)*(ofdm.blockLength + ofdm.CPLength)+(1:(ofdm.blockLength + ofdm.CPLength))) = xcp1;
%     end
% end


%% Plot the PSD of GFDM and OFDM
 f = linspace(-gfdm.K/2, gfdm.K/2, 2*length(sGFDM)+1); 
 f = f(1:end-1)';
 figure(6);
%  plot(f, mag2db(fftshift(abs(fft(sOFDM, 2*length(sOFDM)))))/2, 'b');
 hold on;
 plot(f, mag2db(fftshift(abs(fft(sGFDM, 2*length(sGFDM)))))/2, 'r');
 hold off;
 ylim([-60, 40]);
 xlabel('f/F'); ylabel('PSD [dB]');
 grid()
 legend({'OFDM', 'GFDM'});