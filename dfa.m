function result = dfa(rr, scale_min, scale_max_short, scale_max_long, n_scales)
% DFA  Detrended Fluctuation Analysis for RR interval time series
%
% ECSoC Framework — MATLAB implementation
% Computes α₁ (short-scale), α₂ (long-scale), and R² (goodness-of-fit).
%
% SYNTAX:
%   result = dfa(rr)
%   result = dfa(rr, scale_min, scale_max_short, scale_max_long, n_scales)
%
% INPUTS:
%   rr              — RR interval vector (ms), any orientation
%   scale_min       — Minimum box size (beats). Default: 4
%   scale_max_short — Short-scale upper limit (beats). Default: 16
%   scale_max_long  — Long-scale upper limit (beats). Default: 64
%   n_scales        — Number of log-spaced scales. Default: 20
%
% OUTPUT (struct):
%   .alpha1     — Short-scale DFA exponent
%   .alpha2     — Long-scale DFA exponent
%   .r_squared  — Goodness-of-fit R² (full scale range)
%   .scales     — Box sizes used
%   .F_n        — DFA fluctuation values F(n)
%   .N          — Series length
%
% REFERENCE:
%   Peng CK et al. Chaos 1995;5(1):82–87.
%   Chen Z et al. Phys Rev E 2002;65:041107.

    if nargin < 2, scale_min       = 4;  end
    if nargin < 3, scale_max_short = 16; end
    if nargin < 4, scale_max_long  = 64; end
    if nargin < 5, n_scales        = 20; end

    rr = rr(:);   % force column vector
    N  = numel(rr);

    assert(N >= scale_max_long * 4, ...
        'DFA: series too short. Need >= %d beats, got %d.', scale_max_long * 4, N);
    assert(all(rr > 0), 'DFA: RR intervals must be positive.');

    % Step 1: Cumulative sum (integration)
    y = cumsum(rr - mean(rr));

    % Step 2: Log-spaced box sizes
    raw = unique(round(exp(linspace(log(scale_min), log(scale_max_long), n_scales))));
    scales = raw(raw >= scale_min & raw <= scale_max_long & raw <= floor(N/4));

    assert(numel(scales) >= 4, 'DFA: too few valid scales. Use a longer RR series.');

    % Step 3: Compute F(n) per box size
    F_n = zeros(1, numel(scales));
    for si = 1:numel(scales)
        s       = scales(si);
        n_boxes = floor(N / s);
        sq_res  = zeros(n_boxes, 1);
        for b = 1:n_boxes
            idx  = (b-1)*s + (1:s);
            t    = (1:s)';
            A    = [t, ones(s,1)];
            c    = A \ y(idx);
            sq_res(b) = mean((y(idx) - A*c).^2);
        end
        F_n(si) = sqrt(mean(sq_res));
    end

    % Step 4: Log-log regression per regime
    log_s  = log10(scales);
    log_Fn = log10(F_n);

    idx_s = scales >= scale_min       & scales <= scale_max_short;
    idx_l = scales >= scale_max_short & scales <= scale_max_long;

    p_short = polyfit(log_s(idx_s), log_Fn(idx_s), 1);
    p_long  = polyfit(log_s(idx_l), log_Fn(idx_l), 1);

    % Step 5: Global R²
    p_all   = polyfit(log_s, log_Fn, 1);
    yhat    = polyval(p_all, log_s);
    ss_res  = sum((log_Fn - yhat).^2);
    ss_tot  = sum((log_Fn - mean(log_Fn)).^2);

    result.alpha1    = p_short(1);
    result.alpha2    = p_long(1);
    result.r_squared = 1 - ss_res / ss_tot;
    result.scales    = scales;
    result.F_n       = F_n;
    result.N         = N;
end
