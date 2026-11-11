%% ring_adddrop_fdtd_2d.m
% 2D FDTD (TE: Ex, Ey, Hz) of an ADD-DROP ring resonator:
%   bottom bus waveguide  --gap--  ring  --gap--  top (drop) waveguide
% all embedded in SiO2. Geometry follows the Tidy3D example
% https://www.flexcompute.com/tidy3d/examples/notebooks/RingResonator/
% (R = 5 um, w = 0.43 um, gap = 0.1 um, 1.50 to 1.60 um window).
%
% Light is launched in the bottom bus, travelling +x. The ring is then
% excited counter-clockwise, so the drop signal leaves the top waveguide
% travelling toward -x (left end).
%
% Runs two simulations plus an optional field-map run:
%   REF  : bottom bus only (incident power normalization)
%   RING : bus + ring + drop waveguide
%   MAPS : RING again with full-domain Hz DFT at chosen wavelengths
%
% Outputs: Through = T, Drop = D, Reflection = R (all normalized to the
% incident power), resonance tables (Q, extinction ratio, drop peak).
%f
% Requires MATLAB R2016b or newer. Lengths in micrometers.
% Origin (0,0) is the lower-left corner of the domain.

clear; 
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

%%  USER PARAMETERS -------------------------------
% 2D simulation: use the slab effective index of the 220 nm Si layer
% (about 2.8 to 2.85 for TE at 1.55 um) as n_core, and SiO2 as cladding.
P.n_core = 2.85;
P.n_clad = 1.444;

% Grid and domain
P.dx    = 0.03;            % cell size [um]
P.Lx    = 17;              % domain length in x [um]
P.npml  = 30;              % PML thickness [cells]

% Waveguide and ring geometry
P.wg_w   = 0.43;           % bus / drop waveguide width [um]
P.ring_w = 0.43;           % ring width [um]
P.R      = 5.0;            % ring central radius [um]
P.gap    = 0.10;           % bus-ring and ring-drop gap [um] (set P.gap_drop to differ)
P.gap_drop = P.gap;
P.margin = 1.8;            % clearance between outer waveguide edges and domain edge [um]

% Derived layout (bottom bus, ring above, drop waveguide on top)
P.bus.w  = P.wg_w;
P.bus.cy = P.margin + P.bus.w/2;
P.ring.w  = P.ring_w;
P.ring.R  = P.R;
P.ring.cx = P.Lx/2;
P.ring.cy = P.bus.cy + P.bus.w/2 + P.gap + P.ring.w/2 + P.R;
P.drop.w  = P.wg_w;
P.drop.cy = P.ring.cy + P.R + P.ring.w/2 + P.gap_drop + P.drop.w/2;
P.Ly      = P.drop.cy + P.drop.w/2 + P.margin;

% Source and monitors (x positions in um)
P.x_src  = 1.4;            % soft source plane on the bottom bus
P.x_refl = 2.2;            % reflection monitor (bottom bus, left of ring)
P.x_drop = 2.2;            % drop monitor (top waveguide, left end)
P.x_out  = P.Lx - 2.2;     % through monitor (bottom bus, right of ring)

% Excitation and spectrum
P.lam0      = 1.55;        % center wavelength [um]
P.pulse_tau = 15;          % Gaussian pulse 1/e half-duration [fs]
P.lam_min   = 1.50;        % spectrum window [um]
P.lam_max   = 1.60;
P.nf        = 800;         % number of DFT frequencies
P.T_sim     = 6.0;         % simulated time [ps] (increase for higher Q)
P.dft_stride = 2;          % accumulate DFT every N steps

% Options
P.use_gpu      = false;    % true: use gpuArray (Parallel Computing Toolbox)
P.preview_only = false;    % true: only draw geometry and stop
P.do_maps      = true;     % extra run with full-domain field maps
P.map_lam      = [];       % wavelengths [um] for maps; empty = auto (deepest through dips)
%% END OF INPUTS



%%  GEOMETRY PREVIEW 
Nx = round(P.Lx/P.dx);
Ny = round(P.Ly/P.dx);
xc = ((1:Nx)-0.5)*P.dx;
yc = ((1:Ny)-0.5)*P.dx;
cfg_full = struct('ring', true, 'drop', true);
epsP = eps_map(P, xc, yc, cfg_full);

fprintf('Domain: %.2f x %.2f um, grid %d x %d\n', P.Lx, P.Ly, Nx, Ny);
fprintf('Gap bus-ring  = %.3f um (%.1f cells)\n', P.gap, P.gap/P.dx);
fprintf('Gap ring-drop = %.3f um (%.1f cells)\n', P.gap_drop, P.gap_drop/P.dx);

figure('Name','Geometry');
imagesc(xc, yc, epsP'); axis xy equal tight; colorbar; hold on;
plot([P.x_src P.x_src],[0 P.Ly],'r--');
plot([P.x_refl P.x_refl],[0 P.Ly],'m--');
plot([P.x_out P.x_out],[0 P.Ly],'g--');
dpml = P.npml*P.dx;
plot([dpml P.Lx-dpml P.Lx-dpml dpml dpml],[dpml dpml P.Ly-dpml P.Ly-dpml dpml],'w:');
xlabel('x [\mum]'); ylabel('y [\mum]');
title('\epsilon_r');
drawnow;
print -dpng figure0
if P.preview_only, return; end

%%  SIMULATIONS -----------------------------------
fx = @(E,H) 0.5*real(sum(E.*conj(H),2));      % x-directed Poynting flux per frequency

lam = linspace(P.lam_min, P.lam_max, P.nf).';  % um

fprintf('\n=== REF (bus only) ===\n');
ref = run_fdtd(P, struct('ring',false,'drop',false), [], lam);
fprintf('\n=== ADD-DROP (bus + ring + drop waveguide) ===\n');
r   = run_fdtd(P, cfg_full, [], lam);

A = analyze(ref, r, fx);

%%  RESONANCE TABLES 
tabT = find_resonances(lam, A.T);
tabD = find_drop_peaks(lam, A.D);
print_table_T('Through port dips', tabT);
print_table_D('Drop port peaks', tabD);

%%  PLOTS 
figure('Name','Through and Drop');
subplot(2,1,1);
plot(lam, A.T, 'b', lam, A.D, 'r', lam, A.R, 'k:'); 
grid on;
legend('Through','Drop','Reflection','Location','best');
xlabel('Wavelength [\mum]'); 
ylabel('Power fraction');
title('Add-drop ring: transmission spectra');
%
subplot(2,1,2);
plot(lam, 10*log10(A.T), 'b', lam, 10*log10(A.D), 'r'); 
grid on;
legend('Through','Drop','Location','best');
xlabel('Wavelength [\mum]'); 
ylabel('Transmission [dB]');
ylim([-40 2]);
print -dpng figure1


figure('Name','Power budget');
plot(lam, A.T, 'b', lam, A.D, 'r', lam, A.R, 'k', ...
     lam, 1-A.T-A.D-A.R, 'g'); grid on;
legend('Through','Drop','Reflection','1 - T - D - R (loss, scattering, other ports)','Location','best');
xlabel('Wavelength [\mum]'); ylabel('Power fraction');
title('Power budget');
print -dpng figure2

%%  FIELD MAPS ------------------------------------
if P.do_maps
    if isempty(P.map_lam)
        if isempty(tabT)
            map_lam = P.lam0;
        else
            [~,idx] = sort(tabT(:,2),'ascend');
            map_lam = tabT(idx(1:min(2,numel(idx))),1).';
        end
    else
        map_lam = P.map_lam(:).';
    end
    fprintf('\n=== MAPS at lambda = %s um ===\n', mat2str(map_lam,5));
    rm = run_fdtd(P, cfg_full, map_lam, lam);
    nl = numel(map_lam);
    figure('Name','Field maps');
    for l = 1:nl
        I = abs(rm.Hm(:,:,l)).^2;
        subplot(1,nl,l);
        imagesc(xc, yc, 10*log10(I'/max(I(:)))); axis xy equal tight; hold on;
        contour(xc, yc, epsP', [0.5*(P.n_core^2+P.n_clad^2)], 'w');
        caxis([-40 0]); colorbar;
        xlabel('x [\mum]'); ylabel('y [\mum]');
        title(sprintf('|H_z|^2 [dB], \\lambda = %.4f \\mum', map_lam(l)));
    end
    print -dpng figure3
end

save('ring_adddrop_results.mat','P','lam','A','tabT','tabD');
fprintf('\nDone. Results saved to ring_adddrop_results.mat\n');