function print_table_T(name, tab)
    fprintf('\n%s\n', name);
    if isempty(tab), fprintf('  none found in the window\n'); return; end
    fprintf('  lambda [um]   T_min    ER [dB]   FWHM [nm]   Q_loaded\n');
    for i = 1:size(tab,1)
        fprintf('  %10.4f   %6.3f   %7.2f   %9.3f   %8.0f\n', ...
            tab(i,1), tab(i,2), -10*log10(tab(i,2)), tab(i,3)*1e3, tab(i,4));
    end
end
