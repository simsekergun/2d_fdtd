function tab = find_resonances(lam, T)
    % through-port dips. rows: [lambda_res, T_min, FWHM, loaded Q]
    tab = zeros(0,4);  N = numel(lam);
    for i = 2:N-1
        if T(i) < T(i-1) && T(i) <= T(i+1) && T(i) < 0.9
            lvl = 0.5*(1 + T(i));
            l = i;  while l > 1 && T(l) < lvl, l = l-1; end
            r = i;  while r < N && T(r) < lvl, r = r+1; end
            if l == 1 || r == N, continue; end
            lamL = lam(l)   + (lvl-T(l))  *(lam(l+1)-lam(l))  /(T(l+1)-T(l));
            lamR = lam(r-1) + (lvl-T(r-1))*(lam(r)-lam(r-1))  /(T(r)-T(r-1));
            fw = lamR - lamL;
            tab(end+1,:) = [lam(i), T(i), fw, lam(i)/fw]; 
        end
    end
end
