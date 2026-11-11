function tab = find_drop_peaks(lam, D)
    % drop-port peaks. rows: [lambda_res, D_max, FWHM, loaded Q]
    tab = zeros(0,4);  N = numel(lam);
    for i = 2:N-1
        if D(i) > D(i-1) && D(i) >= D(i+1) && D(i) > 0.05
            lvl = 0.5*D(i);
            l = i;  while l > 1 && D(l) > lvl, l = l-1; end
            r = i;  while r < N && D(r) > lvl, r = r+1; end
            if l == 1 || r == N, continue; end
            lamL = lam(l)   + (lvl-D(l))  *(lam(l+1)-lam(l))  /(D(l+1)-D(l));
            lamR = lam(r-1) + (lvl-D(r-1))*(lam(r)-lam(r-1))  /(D(r)-D(r-1));
            fw = lamR - lamL;
            tab(end+1,:) = [lam(i), D(i), fw, lam(i)/fw]; 
        end
    end
end
