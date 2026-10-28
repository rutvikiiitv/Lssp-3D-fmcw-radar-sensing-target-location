
%% STEP 1: FMCW RADAR PARAMETERS AND CHIRP GENERATION
clear;
clc;
close all;

%% 1. Define radar parameters
c = 3e8;                  % Speed of light (m/s)
fc = 77e9;                % Radar carrier frequency (Hz)
lambda = c/fc;            % Wavelength (m)

B = 1e9;                  % FMCW bandwidth (Hz)
Tc = 50e-6;               % Chirp duration (seconds)
S = B/Tc;                 % Chirp slope (Hz/s)

%% 2. Repeated-chirp parameters for Doppler processing
idleTime = 10e-6;         % Idle time between chirps (seconds)
Tr = Tc + idleTime;       % Chirp repetition interval (seconds)
M = 128;                  % Number of repeated chirps

%% 3. Sampling parameters
% Complex-baseband simulation: we do not directly sample 77 GHz RF.

fs = 5e9;                 % Sampling frequency (Hz)
Ts = 1/fs;                % Sampling interval (seconds)

N = round(Tc*fs);         % Samples per chirp
t = (0:N-1)/fs;           % Time vector

%% 4. Define ideal RF chirp and its baseband equivalent

% Ideal complex RF chirp:
% x_RF(t) = exp(j*(2*pi*fc*t + pi*S*t^2))
%
% We evaluate its phase mathematically. The time vector below is NOT
% a direct 77 GHz RF ADC sampling grid.
phiTxRF = 2*pi*fc*t + pi*S*t.^2;

% Complex-baseband equivalent used by the existing plots and echo model
tx = exp(1j*pi*S*t.^2);

% Instantaneous frequencies
fBaseband = S*t;
fRF = fc + S*t;

%% 5. Display radar parameters
fprintf('====================================\n');
fprintf(' STEP 1: FMCW CHIRP GENERATION\n');
fprintf('====================================\n');

fprintf('Carrier frequency        : %.2f GHz\n',fc/1e9);
fprintf('Wavelength               : %.4f m\n',lambda);
fprintf('Bandwidth                : %.2f GHz\n',B/1e9);
fprintf('Chirp duration           : %.2f us\n',Tc*1e6);
fprintf('Chirp slope              : %.3e Hz/s\n',S);
fprintf('Sampling frequency       : %.2f GHz\n',fs/1e9);
fprintf('Sampling interval        : %.3f ns\n',Ts*1e9);
fprintf('Samples per chirp        : %d\n',N);
fprintf('Chirps per frame          : %d\n',M);
fprintf('Chirp repetition interval: %.2f us\n',Tr*1e6);

%% 6. Plot instantaneous baseband frequency
figure('Name','Step 1 - FMCW Frequency');

plot(t*1e6,fBaseband/1e9,'LineWidth',1.5);
grid on;
xlabel('Time (\mus)');
ylabel('Baseband frequency (GHz)');
title('FMCW Chirp Instantaneous Frequency');

%% 7. Plot real part of the transmitted chirp
figure('Name','Step 1 - Transmitted Chirp');

subplot(2,1,1);
plot(t*1e6,real(tx),'b');
grid on;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Transmitted Chirp - Complex-Baseband Equivalent');

% Zoom into the final 20 ns
zoomDuration = 20e-9;
startTime = Tc - zoomDuration;
idx = find(t >= startTime & t < Tc);

subplot(2,1,2);
plot((t(idx)-startTime)*1e9,real(tx(idx)),'b.-');
grid on;
xlabel('Time from zoom start (ns)');
ylabel('Amplitude');
title('Baseband Chirp - Zoomed Near End');

%% 8. Plot magnitude (optional)
figure('Name','Step 1 - Chirp Magnitude');

plot(t*1e6,abs(tx),'LineWidth',1.5);
grid on;
xlabel('Time (\mus)');
ylabel('Magnitude');
title('Magnitude of Transmitted Chirp');

fprintf('\nStep 1 complete: FMCW chirp generated successfully.\n');
