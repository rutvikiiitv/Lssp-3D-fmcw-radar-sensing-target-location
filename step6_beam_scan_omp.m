%% STEP 6 (REWRITTEN): 4x4 ARRAY RADAR - BEAM SCAN, THEN OMP COMPRESSIVE SENSING
% Pipeline
%   Step 2 scene -> 16-antenna beat data -> range FFT -> Doppler FFT
%   -> BEAM SCAN (coarse theta, phi) -> OMP compressive sensing (refined theta, phi)
%   -> 3D position
%
% Figures (4, plus printed tables)
%   1. 3D scene: RX array + the SAME targets as Step 2
%   2. 3D: actual vs calculated (OMP) target positions
%   3. Beam scan bars (azimuth, elevation) + OMP result, for the targets in targetsToPlot
%   4. Phase difference at every antenna (relative to RX1), one panel per target
%
% Coordinates: x = horizontal, y = forward, z = up
%   theta = azimuth from +y, phi = elevation above the horizontal plane.
% Needs step1_fmcw_chirp.m and step2_random_3d_targets.m in the same folder.
% No toolboxes needed.

clear; clc; close all;


%% 1. TARGETS FROM STEP 2 (same scene is used by every antenna)
run('step2_random_3d_targets.m');
close all;   % remove Step 1/2 figures; delete this line to keep them

needed = {'xTarget','yTarget','zTarget','vx','vy','vz','numTargets', ...
          'reflectivity','c','fc','lambda','S','B','Tc','Tr','M'};
missing = {};
for iv = 1:numel(needed)
    if ~exist(needed{iv},'var'), missing{end+1} = needed{iv}; end %#ok<SAGROW>
end
if ~isempty(missing)
    error('Step 2 / Step 1 did not define: %s',strjoin(missing,', '));
end

% Settings are defined AFTER run(), because Step 1/2 start with 'clear'
% and would wipe anything defined before it.
targetsToPlot = 1;      % which targets get the beam-scan / OMP graphs, e.g. 1:3
arrayDrawScale = 150;   % array is only ~6 mm wide; draw it this many times larger in Figure 1

targetXYZ0     = [xTarget(:), yTarget(:), zTarget(:)];
targetVelocity = [vx(:), vy(:), vz(:)];
nT = numTargets;

%% 2. ARRAY GEOMETRY (x-z plane, half-wavelength spacing)
Nx = 4;  Nz = 4;  nRx = Nx*Nz;
dx = lambda/2;  dz = lambda/2;
xElem = ((0:Nx-1)-(Nx-1)/2)*dx;
zElem = ((0:Nz-1)-(Nz-1)/2)*dz;
rxX = repmat(xElem(:).',Nz,1);  rxX = reshape(rxX.',[],1);   % x varies fastest
rxZ = repmat(zElem(:),1,Nx);    rxZ = reshape(rxZ.',[],1);
arrayY = 0.05;
txPos  = [0 0 0];
rxPos  = [rxX, arrayY*ones(nRx,1), rxZ];

fsBeat   = 10e6;
Nfast    = round(Tc*fsBeat);
tFast    = (0:Nfast-1)'/fsBeat;
slowTime = (0:M-1)*Tr;
snrDb    = 25;

fprintf('===== STEP 6 =====\nTargets: %d | Array: %dx%d | Samples/chirp: %d | Chirps: %d\n', ...
        nT,Nx,Nz,Nfast,M);

%% 3. COMPLEX BEAT SIGNAL AT EVERY ANTENNA
% phi_TX(t) - phi_RX(t) = 2*pi*fc*tau + 2*pi*S*tau*t - pi*S*tau^2, for t >= tau.
% tau uses the exact two-leg path TX -> target -> this RX, which is what
% creates the tiny phase differences between antennas.
beat = complex(zeros(Nfast,M,nRx));
for irx = 1:nRx
    for it = 1:nT
        for m = 1:M
            pos = targetXYZ0(it,:) + targetVelocity(it,:)*slowTime(m);
            tau = (norm(pos-txPos) + norm(pos-rxPos(irx,:)))/c;
            ph  = 2*pi*fc*tau + 2*pi*S*tau*tFast - pi*S*tau^2;
            beat(:,m,irx) = beat(:,m,irx) + reflectivity(it)*(tFast>=tau).*exp(1j*ph);
        end
    end
end
sigPower   = mean(abs(beat(:)).^2);
noisePower = sigPower/10^(snrDb/10);
beat = beat + sqrt(noisePower/2)*(randn(size(beat))+1j*randn(size(beat)));

%% 4. RANGE FFT AND DOPPLER FFT (keep complex values)
wR = reshape(0.5-0.5*cos(2*pi*(0:Nfast-1)'/(Nfast-1)),Nfast,1,1);
RF = fft(beat.*wR,Nfast,1);
nRB = floor(Nfast/2)+1;
RF  = RF(1:nRB,:,:);
rangeAxis = c*((0:nRB-1)'*(fsBeat/Nfast))/(2*S);

wD = reshape(0.5-0.5*cos(2*pi*(0:M-1)/(M-1)),1,M,1);
RD = fftshift(fft(RF.*wD,M,2),2);
velocityAxis = (-floor(M/2):ceil(M/2)-1)/(M*Tr)*lambda/2;
RDmag = abs(RD);

%% 5. FIND EACH TARGET'S RANGE-DOPPLER CELL (validation association)
% Doppler bin first (from known radial velocity), then the strongest range bin
% in that Doppler column. This is a simulation check, not a blind detector.
trueRange = vecnorm(targetXYZ0,2,2);
rb = zeros(nT,1);  db = zeros(nT,1);  estRange = zeros(nT,1);
for it = 1:nT
    vr = targetVelocity(it,:)*targetXYZ0(it,:)'/trueRange(it);
    [~,db(it)] = min(abs(velocityAxis-vr));
    cand = find(rangeAxis>=max(0,trueRange(it)-0.75) & rangeAxis<=trueRange(it)+0.75);
    if isempty(cand)
        [~,rb(it)] = min(abs(rangeAxis-trueRange(it)));
    else
        strength = squeeze(sum(RDmag(cand,db(it),:).^2,3));
        [~,ib] = max(strength);
        rb(it) = cand(ib);
    end
    estRange(it) = rangeAxis(rb(it));
end

%% 6. STAGE 1 - BEAM SCAN (coarse theta, phi over the full field of view)
thetaGrid = -60:1:60;
phiGrid   = -60:1:60;
nTh = numel(thetaGrid);  nPh = numel(phiGrid);
k0  = 2*pi/lambda;
xc  = rxX-mean(rxX);  zc = rxZ-mean(rxZ);

steer = @(thDeg,phDeg) exp(-1j*k0*(xc*(cosd(phDeg(:)).*sind(thDeg(:))).' + ...
                                   zc*sind(phDeg(:)).'))/sqrt(nRx);
[THg,PHg] = meshgrid(thetaGrid,phiGrid);                 % nPh x nTh
Abeam = steer(THg(:),PHg(:));                            % nRx x (nPh*nTh)

trueAz = atan2d(targetXYZ0(:,1),targetXYZ0(:,2));
trueEl = atan2d(targetXYZ0(:,3),hypot(targetXYZ0(:,1),targetXYZ0(:,2)));
bfAz = zeros(nT,1);  bfEl = zeros(nT,1);
beamMap = cell(nT,1);
zCell   = zeros(nRx,nT);

for it = 1:nT
    z = squeeze(RD(rb(it),db(it),:));
    zCell(:,it) = z/max(norm(z),eps);
    s = reshape(abs(Abeam'*zCell(:,it)).^2,nPh,nTh);     % |a' z|^2
    beamMap{it} = s;
    [~,ix] = max(s(:));
    [ip,ia] = ind2sub(size(s),ix);
    bfAz(it) = thetaGrid(ia);   bfEl(it) = phiGrid(ip);
end

%% 7. STAGE 2 - COMPRESSIVE SENSING (OMP on a fine grid around the beam-scan peak)
% Model: z = A*x with x sparse (one non-zero entry per target). The dictionary
% A is built on a 0.25-degree grid within +-12 degrees of the beam-scan result,
% so the beam scan finds the region and OMP sharpens the answer inside it.
winDeg = 12;  fineStep = 0.25;  ompK = 1;
csAz = zeros(nT,1);  csEl = zeros(nT,1);
for it = 1:nT
    thL = max(thetaGrid(1),bfAz(it)-winDeg):fineStep:min(thetaGrid(end),bfAz(it)+winDeg);
    phL = max(phiGrid(1),  bfEl(it)-winDeg):fineStep:min(phiGrid(end),  bfEl(it)+winDeg);
    [THl,PHl] = meshgrid(thL,phL);
    Aloc = steer(THl(:),PHl(:));
    sel  = ompSolve(Aloc,zCell(:,it),ompK,0.05);
    csAz(it) = THl(sel(1));
    csEl(it) = PHl(sel(1));
end

toXYZ = @(R,az,el) [R.*cosd(el).*sind(az), R.*cosd(el).*cosd(az), R.*sind(el)];
bfXYZ = toXYZ(estRange,bfAz,bfEl);
csXYZ = toXYZ(estRange,csAz,csEl);
bfErr = vecnorm(bfXYZ-targetXYZ0,2,2);
csErr = vecnorm(csXYZ-targetXYZ0,2,2);

%% 8. RESULTS TABLE
Target = (1:nT)';
Results = table(Target,trueRange,estRange,trueAz,bfAz,csAz,trueEl,bfEl,csEl,bfErr,csErr, ...
    'VariableNames',{'Target','TrueRange_m','EstRange_m','TrueAz_deg','BeamAz_deg', ...
    'OMP_Az_deg','TrueEl_deg','BeamEl_deg','OMP_El_deg','BeamPosErr_m','OMP_PosErr_m'});
disp(Results);
fprintf('Mean position error: beam scan %.3f m | beam scan + OMP %.3f m\n', ...
        mean(bfErr),mean(csErr));

%% 9. FIGURE 1 - 3D SCENE: RX ARRAY + THE SAME TARGETS AS STEP 2
cols = lines(nT);
f1 = figure('Name','3D scene','Color','w','Position',[60 60 760 620]);
ax = axes(f1); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
set(ax,'Color','w','XColor','k','YColor','k','ZColor','k','GridColor',[0.7 0.7 0.7]);

pad = 0.6*dx;
cx = [-1 1 1 -1]*(max(abs(xElem))+pad)*arrayDrawScale;
cz = [-1 -1 1 1]*(max(abs(zElem))+pad)*arrayDrawScale;
hArr = patch(ax,cx,zeros(1,4),cz,[0.75 0.85 1],'FaceAlpha',0.6,'EdgeColor',[0.2 0.3 0.8]);
scatter3(ax,rxX*arrayDrawScale,zeros(nRx,1),rxZ*arrayDrawScale,22,[0.1 0.2 0.8],'filled');
hTx = scatter3(ax,0,0,0,90,'k','^','filled');
hT  = scatter3(ax,targetXYZ0(:,1),targetXYZ0(:,2),targetXYZ0(:,3),90,cols,'filled');
for it = 1:nT
    plot3(ax,[0 targetXYZ0(it,1)],[0 targetXYZ0(it,2)],[0 targetXYZ0(it,3)], ...
        ':','Color',cols(it,:),'HandleVisibility','off');
    text(ax,targetXYZ0(it,1),targetXYZ0(it,2),targetXYZ0(it,3)+0.4,sprintf('T%d',it), ...
        'Color','k','FontWeight','bold','HorizontalAlignment','center');
end
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(sprintf('Radar and the %d targets from Step 2',nT));
legend([hArr hTx hT],{sprintf('4x4 RX array (drawn %dx larger)',arrayDrawScale), ...
    'TX (array centre)','Targets'},'Location','northeast');
view(40,25);

%% 10. FIGURE 2 - 3D: ACTUAL VS CALCULATED (BEAM SCAN + OMP) POSITIONS
f2 = figure('Name','Actual vs calculated','Color','w','Position',[100 80 760 620]);
ax = axes(f2); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
set(ax,'Color','w','XColor','k','YColor','k','ZColor','k','GridColor',[0.7 0.7 0.7]);
hR = scatter3(ax,0,0,0,100,'k','^','filled');
hA = scatter3(ax,targetXYZ0(:,1),targetXYZ0(:,2),targetXYZ0(:,3),90,cols,'filled');
hC = scatter3(ax,csXYZ(:,1),csXYZ(:,2),csXYZ(:,3),110,cols,'x','LineWidth',2);
for it = 1:nT
    plot3(ax,[targetXYZ0(it,1) csXYZ(it,1)],[targetXYZ0(it,2) csXYZ(it,2)], ...
        [targetXYZ0(it,3) csXYZ(it,3)],'-','Color',cols(it,:),'HandleVisibility','off');
    text(ax,targetXYZ0(it,1),targetXYZ0(it,2),targetXYZ0(it,3)+0.4,sprintf('T%d',it), ...
        'Color','k','FontWeight','bold','HorizontalAlignment','center');
end
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(sprintf('Actual (dot) vs calculated (x) position, mean error %.2f m',mean(csErr)));
legend([hR hA hC],{'Radar','Actual position','Calculated (beam scan + OMP)'}, ...
    'Location','northeast');
view(40,25);

%% 11. FIGURE 3 - BEAM SCAN BARS AND OMP RESULT (one figure per chosen target)
for it = targetsToPlot(targetsToPlot>=1 & targetsToPlot<=nT)
    [~,ipEst] = min(abs(phiGrid-bfEl(it)));
    [~,iaEst] = min(abs(thetaGrid-bfAz(it)));
    azCut = beamMap{it}(ipEst,:);          % azimuth cut at the beam-scan elevation
    elCut = beamMap{it}(:,iaEst).';        % elevation cut at the beam-scan azimuth

    f3 = figure('Name',sprintf('Target %d scan',it),'Color','w','Position',[80 40 1100 720]);
    tl = tiledlayout(f3,2,2,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(tl);
    barScan(ax1,thetaGrid,azCut,bfAz(it),10);
    title(ax1,sprintf('Stage 1: beam scan, azimuth  (true %.1f\\circ)',trueAz(it)));

    ax2 = nexttile(tl);
    barScan(ax2,phiGrid,elCut,bfEl(it),10);
    title(ax2,sprintf('Stage 1: beam scan, elevation  (true %.1f\\circ)',trueEl(it)));

    ax3 = nexttile(tl);
    csPlot(ax3,thetaGrid,azCut,csAz(it),trueAz(it));
    xlabel(ax3,'Azimuth (deg)');
    title(ax3,sprintf('Stage 2: OMP azimuth = %.2f\\circ',csAz(it)));

    ax4 = nexttile(tl);
    csPlot(ax4,phiGrid,elCut,csEl(it),trueEl(it));
    xlabel(ax4,'Elevation (deg)');
    title(ax4,sprintf('Stage 2: OMP elevation = %.2f\\circ',csEl(it)));

    title(tl,sprintf('Target %d: beam scan finds the region, OMP sharpens the angle',it));
end

%% 12. PHASE DIFFERENCE AT EVERY ANTENNA (relative to RX1)
% Measured: phase of the complex range-Doppler value at each antenna minus the
% phase at RX1 (this removes the common ~4*pi*R/lambda offset, leaving only the
% angle-dependent part). Predicted: the same quantity for the OMP angle.
phaseDiffDeg = zeros(nRx,nT);   % measured, wrapped to [-180,180]
predDiffDeg  = zeros(nRx,nT);   % predicted from the OMP (theta,phi)
for it = 1:nT
    zc1 = zCell(:,it);
    phaseDiffDeg(:,it) = rad2deg(angle(zc1.*conj(zc1(1))));
    a = steer(csAz(it),csEl(it));
    predDiffDeg(:,it)  = rad2deg(angle(a.*conj(a(1))));
end
mismatchDeg = rad2deg(angle(exp(1j*deg2rad(phaseDiffDeg-predDiffDeg))));

[ixG,izG] = ndgrid(1:Nx,1:Nz);   % antenna index k -> column ixG(k), row izG(k)
for it = targetsToPlot(targetsToPlot>=1 & targetsToPlot<=nT)
    PhaseTable = table((1:nRx)',ixG(:),izG(:),phaseDiffDeg(:,it),predDiffDeg(:,it), ...
        mismatchDeg(:,it),'VariableNames', ...
        {'Antenna','Column_x','Row_z','Measured_deg','Predicted_OMP_deg','Mismatch_deg'});
    fprintf('\nPhase difference vs RX1, target %d:\n',it);
    disp(PhaseTable);
end

f4 = figure('Name','Phase difference per antenna','Color','w','Position',[100 60 1100 700]);
tl = tiledlayout(f4,'flow','TileSpacing','compact','Padding','compact');
for it = 1:nT
    ax = nexttile(tl);
    P = reshape(phaseDiffDeg(:,it),Nx,Nz).';           % rows = z, columns = x
    imagesc(ax,xElem*1e3,zElem*1e3,P); axis(ax,'xy','equal','tight');
    colormap(ax,hsv); caxis(ax,[-180 180]); hold(ax,'on');
    for k = 1:nRx
        text(ax,xElem(ixG(k))*1e3,zElem(izG(k))*1e3,sprintf('%.0f',phaseDiffDeg(k,it)), ...
            'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color','k');
    end
    set(ax,'XTick',round(xElem*1e3,2),'YTick',round(zElem*1e3,2),'FontSize',7);
    xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
    title(ax,sprintf('T%d  (az %.0f\\circ, el %.0f\\circ)',it,trueAz(it),trueEl(it)));
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'Phase vs RX1 (deg)';
title(tl,'Measured phase difference at each of the 16 antennas (RX1 = bottom-left = 0\circ)');

fprintf('\nStep 6 complete.\n');

%% ===== LOCAL FUNCTIONS =====
function sel = ompSolve(A,z,K,tol)
% Orthogonal Matching Pursuit: choose the atom most correlated with the
% residual, re-fit all chosen atoms by least squares, repeat up to K times.
    r = z;  sel = zeros(1,0);
    for k = 1:K
        [~,j] = max(abs(A'*r));
        sel(end+1) = j;                         %#ok<AGROW>
        coef = A(:,sel)\z;
        r = z - A(:,sel)*coef;
        if norm(r) < tol*norm(z), break; end
    end
end

function barScan(ax,angGrid,cut,bestAng,coarseStep)
% Horizontal bars of normalized score (0-100) at a few candidate angles.
    angs = unique([angGrid(1):coarseStep:angGrid(end), bestAng]);
    vals = interp1(angGrid,cut(:).',angs);
    vals = 100*vals/max(vals);
    n = numel(angs);
    b = barh(ax,1:n,vals,0.55,'FaceColor','flat','EdgeColor','none');
    b.CData = repmat([0.58 0.64 0.74],n,1);
    b.CData(angs==bestAng,:) = [0 0.58 0.5];
    set(ax,'YDir','reverse','YTick',1:n,'YTickLabel',cellstr(compose('%d\\circ',angs)), ...
        'FontSize',8,'Box','off','XGrid','on','GridAlpha',0.15,'Color','w', ...
        'XColor','k','YColor','k');
    xlim(ax,[0 118]); ylim(ax,[0.4 n+0.6]);
    text(ax,vals(:)+2,(1:n)',cellstr(compose('%d',round(vals(:)))),'FontSize',7,'Color','k');
end

function csPlot(ax,angGrid,cut,csAng,trueAng)
% Broad beam-scan curve (grey) vs the single sparse OMP spike (magenta).
    hold(ax,'on'); grid(ax,'on');
    set(ax,'Color','w','XColor','k','YColor','k','Box','off');
    h1 = plot(ax,angGrid,cut(:)/max(cut),'-','Color',[0.6 0.6 0.6],'LineWidth',1.5);
    h2 = stem(ax,csAng,1,'m','filled','LineWidth',2,'MarkerSize',7);
    h3 = plot(ax,[trueAng trueAng],[0 1.1],'k--','LineWidth',1.2);
    ylim(ax,[0 1.15]); ylabel(ax,'Normalized response');
    legend(ax,[h1 h2 h3],{'Beam scan','OMP (sparse)','Truth'},'Location','northeast');
end