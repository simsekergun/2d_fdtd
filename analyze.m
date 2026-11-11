function A = analyze(ref, r, fx)
    Pinc = fx(ref.Ed{2}, ref.Hd{2});                           % incident power (bus only)
    A.T  = fx(r.Ed{2}, r.Hd{2}) ./ Pinc;                       % through port (+x, bottom bus)
    A.D  = -fx(r.Ed{3}, r.Hd{3}) ./ Pinc;                      % drop port (-x, top waveguide)
    A.R  = -fx(r.Ed{1}-ref.Ed{1}, r.Hd{1}-ref.Hd{1}) ./ Pinc;  % reflection (scattered field)
end
