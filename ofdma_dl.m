%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This program illustrates the OFDM signal generation, 
% spectral mask and peak to average power ratio. 
%
%%%%% Parameters %%%%%%%
%
% No.of Carriers: 600 - I.M. for 10 MHz
% Modulation: QPSK, 16QAM or 64 QAM
% Cylic Prefix: Normal
% 
% @author: Ivo Maljevic <ivo.maljevic@utoronto.ca>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ofdma_dl(mod_type)
% mod_type could be 'qpsk', '16qam', '64qam' (you need to be able to run all 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qpsk_table = 1/sqrt(2)*[
    1+j     % 00
    1-j     % 01 
    -1+j    % 10
    -1-j    % 11
]; 

if nargin == 0, mod_type='qpsk'; end

if strcmp(mod_type,'qpsk')
    NBitsPerSymbol = 2;
    Modulator  = comm.PSKModulator(4, 'BitInput', true, 'PhaseOffset', pi/4, 'SymbolMapping', 'Custom', 'CustomSymbolMapping', [0 2 3 1]);
    mod_string = 'QPSK'; 
elseif strcmp(mod_type,'16qam')
    NBitsPerSymbol = 4;
    Modulator = comm.RectangularQAMModulator(16, 'BitInput',true,'NormalizationMethod','Average power', 'SymbolMapping', 'Custom', 'CustomSymbolMapping', ... 
        [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7]);
    mod_string = '16QAM'; 
else
    NBitsPerSymbol = 6;
    Modulator = comm.RectangularQAMModulator(64, 'BitInput',true,'NormalizationMethod','Average power', 'SymbolMapping', 'Custom', 'CustomSymbolMapping', ...
        [47 46 42 43 59 58 62 63 45 44 40 41 ...
        57 56 60 61 37 36 32 33 49 48 52 53 39 ...
        38 34 35 51 50 54 55 7 6 2 3 19 18 22 23 5 ...
        4 0 1 17 16 20 21 13 12 8 9 25 24 28 29 15 ...
        14 10 11 27 26 30 31]);
    mod_string = '64QAM'; 
end

AVG_PWR = 1.3778e-004; % use as is, don't change

N_TTI = 1e4; % specify how many TTIs you want

% End of initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt = 1;
ofdm_data = [];
for n_subframes = 1:N_TTI    % number of subframes
    %fprintf('TTI: %d/%d \n',n_subframes, N_TTI);

    tx_data_symbol = [];
    for n_subcarriers = 1:600
        random_data = randi([0 1], 14*NBitsPerSymbol, 1);           % generate random data for the frame
        modulated_random_data = step(Modulator, random_data);       % modulate random data using QPSK or M-QAM
        tx_data_symbol = [tx_data_symbol; modulated_random_data'];  % order data in a matrix format
    end
    
    for symbol = 1:14                             % cycle through symbols within a TTI
       
        if (symbol == 1 || symbol == 8)           % insert RS pilots and longer CP
            tx_symbol = tx_data_symbol(:,symbol); % one column = one symbol
            % Insert pilots 
            tx_symbol(1:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1));          % set a pilot value here
            tx_symbol(7:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1));                      
	        CP = 159;
        elseif (symbol == 5 || symbol == 12)      % insert RS pilots and normal CP
            tx_symbol = tx_data_symbol(:,symbol); % take one column
            tx_symbol(4:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1));          % RS pilots, need to look up the standard for the more accurate value
            tx_symbol(10:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1));         % RS pilots, need to look up the standard for the more accurate value
            CP = 143; 
        else                                      % normal CP, no RS pilots
            tx_symbol = tx_data_symbol(:,symbol); % take one column
            CP = 143; 
        end
        tx_with_guards = [zeros(723,1); tx_symbol(1:300); 0; tx_symbol(301:end); zeros(724,1)]; % zero DC, other data centered around 0, guards on on each side             
        ofdm_tx_symbol = ifft(tx_with_guards);    % do the IFFT
        ofdm_tx_symbol = [ofdm_tx_symbol(end-CP:end); ofdm_tx_symbol]; % CP insertion at the beginning

        % Peak and average power and their ratio
        peak_power = max(abs(ofdm_tx_symbol).^2); 
        ppr(cnt) = peak_power/AVG_PWR; 
        cnt = cnt +1; 
        
        if n_subframes < 21   % This will speed up the simulation as you don't need all subframes to get the PSD
            ofdm_data = [ofdm_data; ofdm_tx_symbol];   % Just create one large set of data by appending vectors
        end
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PAPR

% Let us plot the CCDF of the peak to average power ratio
bins = sort(ppr);
cdf  = linspace(0,1,length(bins));
ccdf = 1 - cdf;

ylim([1e-2 1])
semilogy(10*log10(bins), ccdf)
grid
t_string = sprintf('PAPR CCDF for %s',mod_string); 
title(t_string);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot PSD

Fs = 30.72e6; % sampling rate
Ts = 1/Fs;    % sampling interval
t_s = 0:Ts:Ts*(length(ofdm_data)-1); 

[Pxx, W] = pwelch(ofdm_data, 2048, 512, [], Fs);  % the trick is how pwelch finds the positive freqs

P_Tx = 46;    % dBm
P_Tx_PSD = P_Tx - 10*log10(9e6); % assuming all power is in 10 MHz
Pxx = Pxx/max(Pxx);   % normalize to 0 dBm

% Spectrum emission mask from the table 6.6.3.1-6 (TS 36.104), f > 1GHz
% case, 5, 10, 15 and 20 MHz channel bandwidths
P_Outside_0 = -7-10*log10(1e5);  % Just the edge
P_Outside_1 = -14-10*log10(1e5); % - 14 dBm / 100 kHz, 5 < df < 10 MHz
P_Outside_2 = -13-10*log10(1e6); % - 13 dBm / 1 MHz, df > 10 MHz

f_sp_mask = [-20,        -15,       -15,         -10,          -5,         -5, 5, 5,           10,         15,           15,           20]*1e6; 
sp_mask = [P_Outside_2, P_Outside_2, P_Outside_1, P_Outside_1,  P_Outside_0, 0, 0, P_Outside_0, P_Outside_1, P_Outside_1, P_Outside_2, P_Outside_2]; 


FS_Half = Fs/2;
Wp = 4.5e6/FS_Half; 
Ws = 5e6/FS_Half;
Rp = 1; 
Rs = 40;

% Calculate the minimum order of Butterworth filter required to meet the
%filter design specifications
[n,Wn] = buttord(Wp,Ws,Rp,Rs); 
% Design the Butterworth filter
[b,a] = butter(n,Wn);

ofdm_data = ofdm_data.*exp(-j*2*pi*Fs/2*t_s');    % shift to the left to prepare for filtering

% Filter the data
ofdm_data_filtered = filter(b,a,ofdm_data).*exp(j*2*pi*Fs/2*t_s'); % filter and shift back for pwelch
[Pxx_filtered, W2] = pwelch(ofdm_data_filtered, 2048, 512, [], Fs);
Pxx_filtered = Pxx_filtered/max(Pxx_filtered);   % normalize to 0 dBm

figure;
plot(W-Fs/2, 10*log10(Pxx)+P_Tx_PSD,f_sp_mask, sp_mask, W-Fs/2, 10*log10(Pxx_filtered)+P_Tx_PSD)
grid
axis auto
t_string = sprintf('Power Spectral Density for %s',mod_string); 
title(t_string);
legend('Signal PSD','SEM','Filtered PSD')
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dBm/Hz)');
