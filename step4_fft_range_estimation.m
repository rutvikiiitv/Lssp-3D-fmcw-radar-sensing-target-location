
%% STEP 4: FFT-BASED RANGE ESTIMATION
% Uses the low-pass-filtered, 10 MHz sampled beat signal from Step 3.
clear;
clc;
close all;

% Run Step 3 to generate the sampled beat signal and target scene
run('step3_beat_signal_sampling.m');

%% 1. Apply a Hann window along fast time
rangeWindow = 0.5 - 0.5*cos( ...
    2*pi*(0:Nfast-1)'/(Nfast-1));

beatWindowed = beat .* rangeWindow;

%% 2. Compute the range FFT
% FFT is calculated along dimension 1 (fast time).
rangeFFT = fft(beatWindowed,Nfast,1);

% Keep positive-frequency bins
numRangeBins = floor(Nfast/2)+1;
rangeFFT = rangeFFT(1:numRangeBins,:);

%% 3. Create frequency and range axes
fRange = (0:numRangeBins-1)'*(fsBeat/Nfast);

% FMCW range conversion
rangeAxis = c*fRange/(2*S);

% Average the magnitude spectrum across all chirps
meanRangeSpectrum = mean(abs(rangeFFT),2);

%% 4. Search for peaks within the expected range
Rmax = 20;

searchBins = find(rangeAxis > 0 & rangeAxis <= Rmax);
searchSpectrum = meanRangeSpectrum(searchBins);

% Find local maxima
peakMask = false(size(searchSpectrum));

if numel(searchSpectrum) >= 3
    peakMask(2:end-1) = ...
        searchSpectrum(2:end-1) > searchSpectrum(1:end-2) & ...
        searchSpectrum(2:end-1) >= searchSpectrum(3:end);
end

candidatePeaks = find(peakMask);

if isempty(candidatePeaks)
    error('No range peaks found. Check the beat signal and FFT.');
end

% Sort candidate peaks from strongest to weakest
[~,order] = sort(searchSpectrum(candidatePeaks),'descend');
candidatePeaks = candidatePeaks(order);

%% 5. Select up to numTargets separated peaks
% The theoretical range resolution is c/(2B).
rangeResolution = c/(2*B);

% Range represented by one FFT bin
rangeBinSpacing = c*(fsBeat/Nfast)/(2*S);

% Approximate minimum separation between selected peak bins
minBins = max(1,round(rangeResolution/rangeBinSpacing));

selectedPeaks = [];

for k = 1:numel(candidatePeaks)

    candidate = candidatePeaks(k);

    if isempty(selectedPeaks) || ...
            all(abs(candidate-selectedPeaks) >= minBins)

        selectedPeaks(end+1) = candidate; %#ok<SAGROW>
    end

    if numel(selectedPeaks) >= numTargets
        break;
    end
end

%% 6. Convert peak frequencies to target ranges
selectedBins = searchBins(selectedPeaks);

estimatedFrequency = fRange(selectedBins);
estimatedRange = c*estimatedFrequency/(2*S);

% Sort estimates from nearest to farthest
[estimatedRange,sortOrder] = sort(estimatedRange);
estimatedFrequency = estimatedFrequency(sortOrder);
selectedBins = selectedBins(sortOrder);

%% 7. Display range estimation results
fprintf('\n====================================\n');
fprintf(' STEP 4: RANGE FFT RESULTS\n');
fprintf('====================================\n');

fprintf('Theoretical range resolution: %.3f m\n',rangeResolution);
fprintf('FFT range-bin spacing: %.3f m\n',rangeBinSpacing);

RangeResults = table( ...
    (1:numel(estimatedRange))', ...
    estimatedFrequency(:), ...
    estimatedRange(:), ...
    'VariableNames', ...
    {'Detection','BeatFrequency_Hz','EstimatedRange_m'});

disp('Detected target candidates:');
disp(RangeResults);

fprintf('Ground-truth ranges from Step 2 (sorted):\n');
disp(sort(R0(:)));

%% 8. Plot the range spectrum
figure('Name','Step 4 - Range FFT Spectrum');

plot(rangeAxis,meanRangeSpectrum,'b','LineWidth',1.2);
hold on;

plot(estimatedRange,meanRangeSpectrum(selectedBins), ...
    'ro','MarkerFaceColor','r');

grid on;
xlim([0 Rmax]);
xlabel('Range (m)');
ylabel('Mean FFT magnitude');
title('Range FFT and Detected Peaks');
legend('Range spectrum','Selected peaks');

%% 9. Compare estimated ranges with true ranges
figure('Name','Step 4 - Range Comparison');

trueRangesSorted = sort(R0(:));

stem(1:numel(trueRangesSorted),trueRangesSorted, ...
    'b','LineWidth',1.5);
hold on;

stem(1:numel(estimatedRange),estimatedRange(:), ...
    'r--','LineWidth',1.5);

grid on;
xlabel('Target index (sorted independently)');
ylabel('Range (m)');
title('True Target Ranges vs FFT Estimates');
legend('True ranges','Estimated ranges');

fprintf('\nStep 4 complete: range estimation finished.\n');
