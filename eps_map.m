function eps = eps_map(P, xv, yv, cfg)
    % relative permittivity sampled at positions xv (1xN) and yv (1xM) in um,
    % with 3x3 supersampling to smooth the staircase of curved interfaces
    ss = 3;  offs = ((1:ss)-(ss+1)/2)/ss*P.dx;
    acc = zeros(numel(xv), numel(yv));
    for a = offs
        for b = offs
            X = xv(:) + a;  Y = yv(:).' + b;
            m = (abs(Y - P.bus.cy) <= P.bus.w/2) & true(numel(xv),1);
            if cfg.drop
                m = m | ((abs(Y - P.drop.cy) <= P.drop.w/2) & true(numel(xv),1));
            end
            if cfg.ring
                rr = sqrt((X - P.ring.cx).^2 + (Y - P.ring.cy).^2);
                m = m | (abs(rr - P.ring.R) <= P.ring.w/2);
            end
            acc = acc + m;
        end
    end
    frac = acc / ss^2;
    eps  = P.n_clad^2 + (P.n_core^2 - P.n_clad^2) * frac;
end
