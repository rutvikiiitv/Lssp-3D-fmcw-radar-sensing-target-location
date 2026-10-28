
%% STEP 5: DOPPLER FFT AND VELOCITY ESTIMATION
% Uses the filtered and sampled beat matrix from Step 3.
% Performs range FFT followed by Doppler FFT.
clear;
clc;
close all;

% Generate the random target scene and repeated-chirp beat signal
run('step3_beat_signal_sampling.m');

%% 1. Apply window along fast time
rangeWindow = 0.5 - 0.5*cos( ...
    2*pi*(0:Nfast-1)'/(Nfast-1));

beatWindowed = beat .* rangeWindow;

%% 2. Perform range FFT
rangeFFT = fft(beatWindowed,Nfast,1);

% Keep positive-frequency range bins
numRangeBins = floor(Nfast/2)+1;
rangeFFT = rangeFFT(1:numRangeBins,:);

fRange = (0:numRangeBins-1)'*(fsBeat/Nfast);
rangeAxis = c*fRange/(2*S);

%% 3. Apply window across repeated chirps
dopplerWindow = 0.5 - 0.5*cos( ...
    2*pi*(0:M-1)/(M-1));

%% 4. Perform Doppler FFT
% Dimension 2 corresponds to the repeated-chirp (slow-time) axis.
rangeDoppler = fftshift( ...
    fft(rangeFFT .* dopplerWindow,M,2),2);

% Doppler frequency axis
fDoppler = (-floor(M/2):ceil(M/2)-1)/(M*Tr);

% Convert Doppler frequency into radial velocity
% Positive velocity means receding in this model.
velocityAxis = fDoppler*lambda/2;

RDmag = abs(rangeDoppler);

%% 5. Plot the range-Doppler map
figure('Name','Step 5 - Range Doppler Map');

imagesc(velocityAxis,rangeAxis, ...
    20*log10(RDmag/max(RDmag(:)) + eps));

axis xy;
ylim([0 20]);

xlabel('Radial velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map');
colorbar;
clim([-40 0]);

%% 6. Estimate radial velocity near each known target range
% This is a validation method: we use ground-truth ranges to select
% nearby range bins. It is not yet a fully automatic target detector.

estimatedVelocity = zeros(numTargets,1);
estimatedRange = zeros(numTargets,1);
estimatedDoppler = zeros(numTargets,1);

for k = 1:numTargets

    % Find the range bin nearest the target's true range
    [~,rangeBin] = min(abs(rangeAxis-R0(k)));

    % Doppler spectrum at this range bin
    dopplerSpectrum = RDmag(rangeBin,:);

    % Select the strongest Doppler bin
    [~,dopplerBin] = max(dopplerSpectrum);

    estimatedDoppler(k) = fDoppler(dopplerBin);
    estimatedVelocity(k) = velocityAxis(dopplerBin);
    estimatedRange(k) = rangeAxis(rangeBin);
end

%% 7. Display velocity results
fprintf('\n====================================\n');
fprintf(' STEP 5: DOPPLER FFT RESULTS\n');
fprintf('====================================\n');

fprintf('Number of chirps: %d\n',M);
fprintf('Chirp repetition interval: %.2f us\n',Tr*1e6);
fprintf('Doppler frequency resolution: %.2f Hz\n',1/(M*Tr));
fprintf('Approximate velocity resolution: %.3f m/s\n', ...
    lambda/(2*M*Tr));
fprintf('Approximate unambiguous velocity: +/- %.3f m/s\n', ...
    lambda/(4*Tr));

DopplerResults = table( ...
    (1:numTargets)', ...
    R0(:), ...
    estimatedRange, ...
    vr(:), ...
    estimatedDoppler, ...
    estimatedVelocity, ...
    isMoving(:), ...
    'VariableNames', ...
    {'Target','TrueRange_m','RangeBin_m', ...
     'TrueRadialVelocity_mps','EstimatedDoppler_Hz', ...
     'EstimatedRadialVelocity_mps','Moving'});

disp(DopplerResults);

%% 8. Plot the Doppler spectrum for one moving target
movingTargets = find(isMoving);

if ~isempty(movingTargets)

    % Choose the moving target with largest absolute radial velocity
    [~,idx] = max(abs(vr(movingTargets)));
    targetToPlot = movingTargets(idx);

    [~,rangeBin] = min(abs(rangeAxis-R0(targetToPlot)));

    figure('Name','Step 5 - Doppler Spectrum');

    plot(velocityAxis,RDmag(rangeBin,:),'LineWidth',1.4);
    hold on;

    grid on;
    xlabel('Radial velocity (m/s)');
    ylabel('FFT magnitude');

    title(sprintf('Doppler Spectrum Near Target %d',targetToPlot));

    xline(vr(targetToPlot),'r--', ...
        'True radial velocity','LineWidth',1.2);
end

fprintf('\nStep 5 complete: Doppler estimation finished.\n');
