
%% STEP 2: RANDOM 3D TARGETS AND RECEIVED ECHOES
clear;
clc;
close all;

% Load radar parameters and transmitted chirp from Step 1
run('step1_fmcw_chirp.m');

%% 1. Generate random 3D target positions
rng('shuffle');

numTargets = 5;

% Radar is located at the origin (0,0,0)
xTarget = -3 + 6*rand(1,numTargets);
yTarget = 3 + 12*rand(1,numTargets);
zTarget = 0.5 + 3*rand(1,numTargets);

% Initial distance from radar to each target
R0 = sqrt(xTarget.^2 + yTarget.^2 + zTarget.^2);

% Target angles relative to radar
azimuth = atan2d(xTarget,yTarget);
elevation = atan2d(zTarget, ...
    sqrt(xTarget.^2 + yTarget.^2));

%% 2. Assign random target motion
% Ensure at least one moving and one stationary target
numMoving = randi([1,numTargets-1]);

isMoving = false(1,numTargets);
movingOrder = randperm(numTargets,numMoving);
isMoving(movingOrder) = true;

% Random velocity components in m/s
vx = -3 + 6*rand(1,numTargets);
vy = -3 + 6*rand(1,numTargets);
vz = -1 + 2*rand(1,numTargets);

% Stationary targets have zero velocity
vx(~isMoving) = 0;
vy(~isMoving) = 0;
vz(~isMoving) = 0;

% Radial velocity: positive = moving away
vr = (xTarget.*vx + yTarget.*vy + zTarget.*vz)./R0;

%% 3. Assign target reflectivity and phase
targetAmplitude = 0.3 + 0.7*rand(1,numTargets);
targetPhase = 2*pi*rand(1,numTargets);

reflectivity = targetAmplitude .* exp(1j*targetPhase);

%% 4. Calculate round-trip delays

% Round-trip propagation delay for each target
tau = 2*R0/c;

%% 5. Generate received baseband echoes for one chirp

% Simplified complex-baseband echo illustration.
% The carrier-dependent propagation phase is included, but the RF
% carrier is not directly sampled.

rx = complex(zeros(size(tx)));

for k = 1:numTargets

    % Receiver time corresponding to the delayed transmitted signal
    delayedTime = t - tau(k);

    % Only add the echo where the delayed chirp exists
    valid = delayedTime >= 0 & delayedTime < Tc;

    % Carrier-dependent propagation phase
    carrierPhase = exp(-1j*2*pi*fc*tau(k));

    % Delayed complex-baseband chirp
    delayedChirp = exp(1j*pi*S*delayedTime(valid).^2);

    % Add this target's echo
    rx(valid) = rx(valid) + ...
        reflectivity(k)*carrierPhase*delayedChirp;

end

%% 6. Display target information
TargetTable = table( ...
    (1:numTargets)', ...
    xTarget',yTarget',zTarget', ...
    R0',azimuth',elevation', ...
    isMoving',vx',vy',vz',vr', ...
    'VariableNames', ...
    {'Target','X_m','Y_m','Z_m','Range_m', ...
     'Azimuth_deg','Elevation_deg','Moving', ...
     'Vx_mps','Vy_mps','Vz_mps','RadialVelocity_mps'});

fprintf('\n====================================\n');
fprintf(' STEP 2: RANDOM 3D TARGETS\n');
fprintf('====================================\n');

disp(TargetTable);

fprintf('Number of targets: %d\n',numTargets);
fprintf('Moving targets: %d\n',sum(isMoving));
fprintf('Stationary targets: %d\n',sum(~isMoving));

%% 7. Plot the 3D target scene
figure('Name','Step 2 - Random 3D Targets');

scatter3(xTarget(~isMoving), ...
    yTarget(~isMoving),zTarget(~isMoving), ...
    100,'b','filled');

hold on;

scatter3(xTarget(isMoving), ...
    yTarget(isMoving),zTarget(isMoving), ...
    100,'r','filled');

scatter3(0,0,0,160,'k','^','filled');

for k = 1:numTargets
    text(xTarget(k),yTarget(k),zTarget(k), ...
        sprintf(' T%d',k),'FontSize',10);
end

grid on;
axis equal;

xlabel('X position (m)');
ylabel('Y position (m)');
zlabel('Z position (m)');

title('Random 3D Radar Target Scene');

legend('Stationary targets','Moving targets', ...
    'Radar','Location','best');

view(40,25);

%% 8. Plot transmitted and received echo magnitudes
figure('Name','Step 2 - Transmitted and Received Signals');

subplot(2,1,1);
plot(t*1e6,abs(tx),'b');
grid on;
xlabel('Time (\mus)');
ylabel('Magnitude');
title('Transmitted Chirp Magnitude');

subplot(2,1,2);
plot(t*1e6,abs(rx),'r');
grid on;
xlabel('Time (\mus)');
ylabel('Magnitude');
title('Combined Received Echo Magnitude');

%% 9. Plot real parts over a short time window
figure('Name','Step 2 - Waveform Comparison');

numShow = min(round(0.1e-6*fs),N);

plot(t(1:numShow)*1e9,real(tx(1:numShow)),'b');
hold on;
plot(t(1:numShow)*1e9,real(rx(1:numShow)),'r');

grid on;
xlabel('Time (ns)');
ylabel('Real amplitude');
title('Transmitted Chirp and Combined Echoes');
legend('Transmitted signal','Received echoes');

fprintf('\nStep 2 complete: random targets and echoes generated.\n');
