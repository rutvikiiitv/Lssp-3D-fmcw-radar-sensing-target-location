
%% STEP 3: RF MIXER, LOW-PASS FILTER, AND BEAT SIGNAL SAMPLING
clear;
clc;
close all;

% Run Step 2 to generate a new random target scene
run('step2_random_3d_targets.m');

%% 1. Define mixer, filter, and ADC sampling parameters

% Final ADC sampling rate for the low-frequency beat signal
fsBeat = 10e6;
TsBeat = 1/fsBeat;

Nfast = round(Tc*fsBeat);
tFast = (0:Nfast-1)'/fsBeat;

% Internal simulation rate for the low-frequency mixer output.
% This is NOT the 77 GHz RF sampling rate.
fsMix = 100e6;

% Filter at fsMix, then downsample to fsBeat
decimationFactor = round(fsMix/fsBeat);

if abs(fsMix/decimationFactor - fsBeat) > 1
    error('Internal rate must be an integer multiple of ADC rate.');
end

Nmix = round(Tc*fsMix);
tMix = (0:Nmix-1)'/fsMix;

% One slow-time instant per chirp
slowTime = (0:M-1)*Tr;

%% 2. Ideal RF phase model and complex mixer

% Ideal transmitted RF phase:
% phi_TX(t) = 2*pi*fc*t + pi*S*t^2
%
% Delayed received RF phase:
% phi_RX(t) = 2*pi*fc*(t-tau) + pi*S*(t-tau)^2
%
% Ideal complex mixing calculates:
% TX * conj(RX) = exp(j*(phi_TX - phi_RX))
%
% We calculate the phase difference directly.
% We do NOT directly sample the 77 GHz carrier.
%
% The changing target range across chirps represents target motion.

mixed = complex(zeros(Nmix,M));

for k = 1:numTargets

    for m = 1:M

        % Target range during this chirp
        Rm = R0(k) + vr(k)*slowTime(m);

        % Round-trip propagation delay
        tau_m = 2*Rm/c;

        % Ideal transmitted RF phase
        phiTx = 2*pi*fc*tMix + pi*S*tMix.^2;

        % Delayed echo time
        delayedTime = tMix - tau_m;

        % Echo exists only after the round-trip delay
        valid = delayedTime >= 0 & delayedTime < Tc;

        % Delayed received RF phase
        phiRx = 2*pi*fc*delayedTime(valid) + ...
                pi*S*delayedTime(valid).^2;

        % Difference phase from ideal complex mixing
        phaseDifference = phiTx(valid) - phiRx;

        % Create this target's mixer output
        targetMixed = complex(zeros(Nmix,1));

        targetMixed(valid) = reflectivity(k) .* ...
            exp(1j*phaseDifference);

        % Add the target contribution
        mixed(:,m) = mixed(:,m) + targetMixed;

    end
end

%% 3. Low-pass filter the mixer output

% Maximum beat frequency for the range of interest
Rmax = 20;
fBeatMax = 2*S*Rmax/c;

% LPF cutoff frequency
lpfCutoff = 3.5e6;

if lpfCutoff <= fBeatMax
    error('LPF cutoff must exceed the maximum expected beat frequency.');
end

filteredMixed = complex(zeros(size(mixed)));

% Toolbox-free first-order digital low-pass filter.
%
% Difference equation:
% y[n] = alpha*y[n-1] + (1-alpha)*x[n]
%
% alpha = exp(-2*pi*cutoff_frequency/sampling_frequency)
%
% This is a first-order RC-style digital LPF approximation.

alphaLPF = exp(-2*pi*lpfCutoff/fsMix);

for m = 1:M

    x = mixed(:,m);
    y = complex(zeros(size(x)));

    % Initialize filter output
    y(1) = (1-alphaLPF)*x(1);

    % Apply the low-pass filter
    for n = 2:Nmix

        y(n) = alphaLPF*y(n-1) + ...
               (1-alphaLPF)*x(n);

    end

    filteredMixed(:,m) = y;

end

%% 4. ADC sampling at 10 MHz

% Downsample only AFTER low-pass filtering
beat = filteredMixed(1:decimationFactor:end,:);

% Ensure exactly Nfast samples per chirp
beat = beat(1:Nfast,:);

%% 5. Add complex Gaussian noise

snrDb = 25;

signalPower = mean(abs(beat(:)).^2);
noisePower = signalPower/10^(snrDb/10);

noise = sqrt(noisePower/2) * ...
    (randn(size(beat)) + 1j*randn(size(beat)));

beat = beat + noise;

%% 6. Sampling Theorem calculations

% Maximum beat frequency for the range of interest
nyquistBeat = 2*fBeatMax;

fprintf('\n====================================\n');
fprintf(' STEP 3: MIXER, LPF, AND ADC SAMPLING\n');
fprintf('====================================\n');

fprintf('RF carrier frequency: %.2f GHz\n',fc/1e9);
fprintf('Internal mixer simulation rate: %.2f MHz\n',fsMix/1e6);
fprintf('LPF cutoff frequency: %.2f MHz\n',lpfCutoff/1e6);
fprintf('Beat ADC sampling frequency: %.2f MHz\n',fsBeat/1e6);
fprintf('ADC sampling interval: %.2f ns\n',TsBeat*1e9);
fprintf('Samples per chirp: %d\n',Nfast);
fprintf('Number of chirps: %d\n',M);

fprintf('Maximum beat frequency at %.1f m: %.3f MHz\n', ...
    Rmax,fBeatMax/1e6);

fprintf('Nyquist rate for beat signal: %.3f MHz\n', ...
    nyquistBeat/1e6);

if fsBeat > nyquistBeat
    fprintf('Sampling check: PASS\n');
else
    fprintf('Sampling check: WARNING - rate may be insufficient\n');
end

%% 7. Create Sampling Theorem table

Signal = [
    "77 GHz RF carrier"
    "Complex-baseband chirp"
    "Beat signal up to 20 m"
];

MaximumFrequency_Hz = [
    fc
    B
    fBeatMax
];

NyquistRate_Hz = 2*MaximumFrequency_Hz;

ChosenSamplingRate_Hz = [
    NaN
    NaN
    fsBeat
];

Comment = [
    "RF phase model; not sampled directly"
    "Equivalent baseband representation"
    "10 MHz ADC samples filtered beat signal"
];

SamplingTable = table(Signal,MaximumFrequency_Hz, ...
    NyquistRate_Hz,ChosenSamplingRate_Hz,Comment);

disp(' ');
disp('===== SAMPLING THEOREM TABLE =====');
disp(SamplingTable);

%% 8. Plot sampled beat signal from the first chirp

figure('Name','Step 3 - Sampled Beat Signal');

plot(tFast*1e6,real(beat(:,1)),'b');
grid on;

xlabel('Fast time (\mus)');
ylabel('Real amplitude');
title('Low-Pass-Filtered Beat Signal: First Chirp');

%% 9. Plot beat-signal magnitude

figure('Name','Step 3 - Beat Magnitude');

plot(tFast*1e6,abs(beat(:,1)),'r');
grid on;

xlabel('Fast time (\mus)');
ylabel('Magnitude');
title('Beat Signal Magnitude: First Chirp');

%% 10. Show individual sample points

figure('Name','Step 3 - Sampling Points');

numShow = min(150,Nfast);

plot(tFast(1:numShow)*1e6, ...
    real(beat(1:numShow,1)),'b-');

hold on;

stem(tFast(1:numShow)*1e6, ...
    real(beat(1:numShow,1)),'r.', ...
    'MarkerSize',8);

grid on;

xlabel('Fast time (\mus)');
ylabel('Real amplitude');
title('Beat Signal and Discrete Sample Values');
legend('Sampled signal','Sample values');

%% 11. Number of samples and sampling information

samplesPerChirp = Nfast;
totalSamples = Nfast*M;
samplingInterval = 1/fsBeat;

fprintf('\n===== SAMPLE COUNT =====\n');
fprintf('Sampling rate: %.3f MHz\n',fsBeat/1e6);
fprintf('Sampling interval: %.3f microseconds\n', ...
    samplingInterval*1e6);
fprintf('Samples per chirp: %d\n',samplesPerChirp);
fprintf('Number of chirps: %d\n',M);
fprintf('Total beat samples: %d\n',totalSamples);

%% 12. Sampling rate diagram

% Display the first 2 microseconds of the first chirp
figure('Name','Step 3 - Sampling Rate Diagram');

displayDuration = 2e-6;
numDisplay = min(round(displayDuration*fsBeat),Nfast);

timeDisplay = tFast(1:numDisplay)*1e6;
signalDisplay = real(beat(1:numDisplay,1));

plot(timeDisplay,signalDisplay,'b-o', ...
    'LineWidth',1.1,'MarkerSize',3);

hold on;

stem(timeDisplay,signalDisplay,'r.', ...
    'LineWidth',1,'MarkerSize',10);

grid on;

xlabel('Time (\mus)');
ylabel('Amplitude');
title('Sampling Rate Diagram of the Beat Signal');
legend('Sample values','Sampling instants');

%% 13. Sampling theorem rate comparison

nyquistRate = 2*fBeatMax;

figure('Name','Step 3 - Sampling Rate Comparison');

bar([nyquistRate,fsBeat]/1e6);

set(gca,'XTickLabel', ...
    {'Nyquist rate','Chosen sampling rate'});

ylabel('Sampling frequency (MHz)');
title('Sampling Theorem: Required vs Chosen Rate');
grid on;

fprintf('\nStep 3 sampling analysis complete.\n');
