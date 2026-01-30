clc; clear; close all;
addpath(genpath('C:\Users\mittaus\Desktop\Dedar'));
addpath(genpath('C:\Users\mittaus\Desktop\Dedar\iqtools'))
%% ======= Connection (edit if needed) =======
Host    = 'localhost';      % if MATLAB runs on the AXIe controller; else use controller host/IP, e.g. 'DESKTOP-69KH8GQ'
LANName = 'hislip0';              % from AgM8190Firmware (e.g., hislip2)
visaRsrc = sprintf('TCPIP0::%s::%s::INSTR', Host, LANName);
psg_device = 0;
uxr_device = 0;
vsa_soft = 1;
connect_2_PSG_UXA_VSA;
fopen(vsa)

%% ======= 64-QAM parameters =======
M           = 64;                 % Modulation order
symR        = 100e6;              % Symbol rate [symbols/s]
frameLength = 16384;             % Symbols per frame
nSamps      = 5;                 % Samples per symbol (baseband)
Fs_bb       = symR * nSamps;      % Baseband sample rate
beta        = 0.35;               % RRC rolloff
filtergain  = 0.4;                % RRC gain
nsymb       = 32;                 % RRC span [symbols]

%% ======= Build baseband waveform =======
data = randi([0 M-1], frameLength, 1);
tx   = qammod(data, M, 'UnitAveragePower', true);

rcFilt = comm.RaisedCosineTransmitFilter( ...
    'RolloffFactor',          beta, ...
    'FilterSpanInSymbols',    nsymb, ...
    'OutputSamplesPerSymbol', nSamps, ...
    'Gain',                   filtergain);

tx_upd    = [tx; zeros(nsymb,1)];
filt_data = step(rcFilt, tx_upd);
Inp_Sig   = filt_data(nsymb*(nSamps/2)+1 : frameLength*nSamps + nsymb*(nSamps/2), 1);

% ===== Symbol alignment (prevents EVM blow-up) =====
Lsym = floor(length(Inp_Sig)/nSamps) * nSamps;
Inp_Sig = Inp_Sig(1:Lsym);

% ===== AWG granularity requirement =====
gran_8190A = 48;      % for M8190A
remg = mod(length(Inp_Sig), gran_8190A);
if remg ~= 0
    Inp_Sig = [Inp_Sig; zeros(gran_8190A - remg, 1)];
end


Inp_QAM_signal = Inp_Sig;
%% ======= Resample to 8 GSa/s & upconvert to 2 GHz =======
Fs_arb1 = 8e9;                     % DAC sample rate (12G option present)
Fs     = Fs_arb1;                  % iqdownload expects "Fs"
arb1_in_at_bb = resample(Inp_Sig, Fs_arb1, Fs_bb);
% arb1_in_at_8gbb = resample(Inp_Sig, Fs_bb);

fIF = 2.4e9;                        % desired center frequency
n   = length(arb1_in_at_bb);
t   = (0:n-1).' / Fs;

arb1_in = arb1_in_at_bb .* exp(1j*2*pi*fIF*t);   % numerical upconversion

n   = length(arb1_in);
%% ======= Marker & channel map =======
marker = [zeros(floor(n/2),1); 15*ones(n-floor(n/2),1)];


channelmapping = [1 0; 0 0];                               % CH1 & CH2


%% ======= IQTools ARB configuration (single M8190A) =======
arbConfig                    = struct;
arbConfig.model              = 'M8190A_14bit';
arbConfig.connectionType     = 'visa';
arbConfig.visaAddr           = visaRsrc;        % TCPIP0::<host>::hislipX::INSTR
arbConfig.visaVendor         = 'keysight';      % force Keysight VISA
arbConfig.defaultFc          = 0;
arbConfig.tooltips           = 1;
arbConfig.DACRange           = 1;
arbConfig.amplScale          = 2.83;
arbConfig.amplScaleMode      = 'Leave Unchanged';
arbConfig.triggerMode        = 'Continuous';
arbConfig.useM8192A          = 0;
arbConfig.isScopeConnected   = 0;
arbConfig.isVSAConnected     = 0;
arbConfig.isDCAConnected     = 0;
arbConfig.OutputBufferSize   = 1e8;
arbConfig.TimeOut            = 20;


arb1 = instrfind('Type', 'visa-tcpip', 'RsrcName', 'TCPIP0::localhost::5025::SOCKET', 'Tag', '');
if isempty(arb1)
    arb1 = tcpip('localhost', 5025);
else
    fclose(arb1);
    arb1 = arb1(1);
end

set(arb1, 'TimeOut', 100);
fopen(arb1);
idn_arb1 = query(arb1, '*IDN?');
fprintf(arb1, ':ROSC:SOUR EXT');
fprintf(arb1, ':ROSC:FREQ 1E7');
% fprintf(arb1, ':MARK:SAMP:VOLT:LOW 0.1');
fprintf(1, 'Matlab connected to %s\n', idn_arb1);
fclose(arb1);

%% ======= Download to M8190A CH1 & CH2 (IQTools, DC arg-free) =======
iqdownload(arb1_in, Fs, ...
    'arbConfig',      arbConfig, ...
    'segmentNumber',  1, ...
    'channelmapping', channelmapping, ...
    'marker',         marker, ...
    'keepOpen',       false);

fprintf('M8190A CH1/CH2 loaded @ %.1f GSa/s | Center = %.1f GHz | Length = %d (×48)\n', Fs/1e9, fIF/1e9, length(arb1_in));


%% VSA Data 
tracenum = 3;
fprintf(vsa, sprintf(':INP:ANAL:RANG:AUTO'));
pause(5);
fprintf(vsa, sprintf(':INIT:REST'));
pause(0.5);
fprintf(vsa,':INIT:IMM;OPC?');
pause(5);
fclose(vsa);
Ns = numel(Inp_QAM_signal);
[iqdata_vsa1] = My_testVSAdata_shun(vsa, Ns,tracenum); % taking data out of VSA

%% Time synchronization of Tx data (inp_signal) and Rx data (iqdata_vsa)

% Find the delay between the two signals using finddelay

% Find the delay between the input signal and the first output signal
a1 = finddelay(Inp_QAM_signal, iqdata_vsa1);
timeDelay1 = 1/Fs * a1;
fprintf('Synchronized time delay of %d s which corresponds to %i samples\n', timeDelay1, a1);

if a1 >= 0
    iqdata_vsa1 = iqdata_vsa1((a1+1):end);
elseif a1<0
    Inp_QAM_signal = Inp_QAM_signal((-a1+1):end);
end

% % % Truncation of Rx data to match with Tx data in size
if numel(iqdata_vsa1) >= numel(Inp_QAM_signal)
    iqdata_vsa1 = iqdata_vsa1(1:numel(Inp_QAM_signal));
end

% Tx and Rx data after time synchronization
xin1 = Inp_QAM_signal(1:numel(iqdata_vsa1)); % iqdata1(1:numel(iqdata1));
Outputsignal = iqdata_vsa1;


% Verify signal synchronization
synch1 = ~finddelay(Inp_QAM_signal, Outputsignal);
xophs1 = xcorr(Inp_QAM_signal, Outputsignal);


figure(2);
plot(abs(xophs1));
title('Cross-correlation of synchronized signals (iqdata\_vsa1)');

%save('C:\Users\mittaus\Desktop\Dedar\output_waveform\Inp_QAM_signal.mat')
