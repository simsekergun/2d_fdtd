function out = run_fdtd(P, cfg, map_lam, lam)
    c0 = 299792458; mu0 = 4e-7*pi; eps0 = 1/(mu0*c0^2);
    dx = P.dx*1e-6;
    Nx = round(P.Lx/P.dx);
    Ny = round(P.Ly/P.dx);
    % stability: limited by the fastest medium (cladding)
    dt = 0.99*dx*P.n_clad/(c0*sqrt(2));
    Nt = ceil(P.T_sim*1e-12/dt);
    G = @(a) a;  if P.use_gpu, G = @gpuArray; end

    % Yee grid (TE): Hz at cell centers, Ex at (center x, node y), Ey at (node x, center y)
    xc = ((1:Nx)-0.5)*P.dx;
    yc = ((1:Ny)-0.5)*P.dx;
    xn = (0:Nx)*P.dx;
    yn = (0:Ny)*P.dx;
    epsEx = eps_map(P, xc, yn, cfg);     % Nx x (Ny+1)
    epsEy = eps_map(P, xn, yc, cfg);     % (Nx+1) x Ny

    % Split-field PML with loss rate kappa [1/s]: sigma = kappa*eps, sigma* = kappa*mu0
    d = P.npml*dx;
    mpow = 3;
    R0 = 1e-6;
    kmax = -(mpow+1)*log(R0)*(c0/P.n_clad)/(2*d);
    Lxm = Nx*dx;  Lym = Ny*dx;
    kap = @(pos,L) kmax*min(1, max(0, max(d-pos, pos-(L-d))/d)).^mpow;
    kx_c = kap(xc*1e-6, Lxm);  kx_c = kx_c(:);       % Nx x 1
    kx_n = kap(xn*1e-6, Lxm);  kx_n = kx_n(:);       % (Nx+1) x 1
    ky_c = kap(yc*1e-6, Lym);  ky_c = ky_c(:).';     % 1 x Ny
    ky_n = kap(yn*1e-6, Lym);  ky_n = ky_n(:).';     % 1 x (Ny+1)
    q = @(k) k*dt/2;

    cax = G(((1-q(ky_n))./(1+q(ky_n))) .* ones(Nx,1));
    cbx = G((dt/(eps0*dx)) ./ (epsEx .* (1+q(ky_n))));
    cay = G(((1-q(kx_n))./(1+q(kx_n))) .* ones(1,Ny));
    cby = G((dt/(eps0*dx)) ./ (epsEy .* (1+q(kx_n))));
    hax = G(((1-q(kx_c))./(1+q(kx_c))) .* ones(1,Ny));
    hbx = G(((dt/(mu0*dx)) ./ (1+q(kx_c))) .* ones(1,Ny));
    hay = G(((1-q(ky_c))./(1+q(ky_c))) .* ones(Nx,1));
    hby = G(((dt/(mu0*dx)) ./ (1+q(ky_c))) .* ones(Nx,1));
    clear epsEx epsEy

    Ex  = G(zeros(Nx,Ny+1));
    Ey  = G(zeros(Nx+1,Ny));
    Hzx = G(zeros(Nx,Ny));
    Hzy = G(zeros(Nx,Ny));
    Hz  = G(zeros(Nx,Ny));

    % Source: soft Hz line source with a Gaussian transverse profile on the bottom bus
    f0 = c0/(P.lam0*1e-6);  
    w0 = 2*pi*f0;
    tau = P.pulse_tau*1e-15;
    t0 = 5*tau;
    gsrc = @(t) exp(-((t-t0)/tau)^2)*cos(w0*(t-t0));
    isrc = round(P.x_src/P.dx + 0.5);
    jsrc = round((P.bus.cy-P.bus.w/2-1.0)/P.dx+0.5):round((P.bus.cy+P.bus.w/2+1.0)/P.dx+0.5);
    prof = G(exp(-((yc(jsrc)-P.bus.cy)/(0.6*P.bus.w)).^2));

    % Monitors (vertical lines). 1: reflection, 2: through, 3: drop
    jbus  = round((P.bus.cy -P.bus.w/2 -0.8)/P.dx+0.5):round((P.bus.cy +P.bus.w/2 +0.8)/P.dx+0.5);
    jdrop = round((P.drop.cy-P.drop.w/2-0.8)/P.dx+0.5):round((P.drop.cy+P.drop.w/2+0.8)/P.dx+0.5);
    mon(1).i = round(P.x_refl/P.dx)+1;  mon(1).jr = jbus;
    mon(2).i = round(P.x_out /P.dx)+1;  mon(2).jr = jbus;
    mon(3).i = round(P.x_drop/P.dx)+1;  mon(3).jr = jdrop;
    nm = 3;

    wv = 2*pi*c0./(lam(:)*1e-6);                    % Nf x 1 angular frequencies
    Ed = cell(1,nm); Hd = cell(1,nm);
    for m = 1:nm
        Ed{m} = G(complex(zeros(numel(wv), numel(mon(m).jr))));
        Hd{m} = Ed{m};
    end

    domaps = ~isempty(map_lam);
    if domaps
        wm = 2*pi*c0./(map_lam(:)*1e-6);
        nl = numel(wm);
        Hm = G(complex(zeros(Nx,Ny,nl)));
    else
        Hm = [];
    end

    fprintf('Grid %d x %d, dt = %.3g fs, %d steps\n', Nx, Ny, dt*1e15, Nt);
    tic;  
    hpeak = 0;  
    hlast = 0;
    for n = 1:Nt
        % E update (uses Hz at n-1/2)
        Ex(:,2:Ny)  = cax(:,2:Ny).*Ex(:,2:Ny) + cbx(:,2:Ny).*(Hz(:,2:Ny)-Hz(:,1:Ny-1));
        Ey(2:Nx,:)  = cay(2:Nx,:).*Ey(2:Nx,:) - cby(2:Nx,:).*(Hz(2:Nx,:)-Hz(1:Nx-1,:));
        % H update (gives Hz at n+1/2)
        Hzx = hax.*Hzx - hbx.*(Ey(2:Nx+1,:)-Ey(1:Nx,:));
        Hzy = hay.*Hzy + hby.*(Ex(:,2:Ny+1)-Ex(:,1:Ny));
        Hzx(isrc,jsrc) = Hzx(isrc,jsrc) + prof*gsrc((n+0.5)*dt);
        Hz = Hzx + Hzy;

        % DFT of monitor lines
        if mod(n,P.dft_stride) == 0
            eE = G(exp(-1i*wv*(n*dt)));
            eH = G(exp(-1i*wv*((n+0.5)*dt)));
            for m = 1:nm
                i = mon(m).i;
                jr = mon(m).jr;
                eyl = Ey(i,jr);
                hzl = 0.5*(Hz(i-1,jr) + Hz(i,jr));
                Ed{m} = Ed{m} + eE*eyl;
                Hd{m} = Hd{m} + eH*hzl;
            end
        end
        % full-domain Hz DFT for field maps
        if domaps && mod(n,4) == 0
            for l = 1:nl
                Hm(:,:,l) = Hm(:,:,l) + exp(-1i*wm(l)*(n+0.5)*dt)*Hz;
            end
        end
        if mod(n,250) == 0
            hlast = gather(max(abs(Hz(:))));  hpeak = max(hpeak, hlast);
        end
        if mod(n,5000) == 0
            el = toc;
            fprintf('  step %d/%d, elapsed %.0f s, ETA %.0f s\n', n, Nt, el, el*(Nt-n)/n);
        end
    end
    fprintf('Finished in %.0f s. Residual field / peak field = %.2e\n', toc, hlast/hpeak);
    if hlast/hpeak > 1e-3
        warning('Fields have not decayed: increase P.T_sim for accurate spectra (high-Q rings).');
    end

    out.Ed = cellfun(@gather, Ed, 'UniformOutput', false);
    out.Hd = cellfun(@gather, Hd, 'UniformOutput', false);
    out.Hm = gather(Hm);
    out.dt = dt;
end