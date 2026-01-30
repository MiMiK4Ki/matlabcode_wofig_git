function [h_butterworth, butterworth_t] = generate_butterworth_filter(f_cutoff, hFLpS, Tc, deltat)
    [b, a] = butter(6, 2 * pi * f_cutoff, 's');
    [resi, p] = residue(b, a);
    butterworth_t = 0:deltat:10 * hFLpS * Tc;
    h_butterworth = real(resi.' * exp(p .* butterworth_t));
end
