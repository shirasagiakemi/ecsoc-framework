function result = run_ecsoc(rr_ms, r2_threshold)
% RUN_ECSOC  Full ECSoC analysis pipeline for a single RR interval series
%
% ECSoC Framework — MATLAB implementation
% Calls: dfa.m → chi_calculation.m → phase_classification.m
%
% SYNTAX:
%   result = run_ecsoc(rr_ms)
%   result = run_ecsoc(rr_ms, r2_threshold)
%
% QUICK START:
%   data   = readmatrix('data/sample/sample_rr_normal.csv');
%   rr     = data(:, 1);   % rr_ms column
%   result = run_ecsoc(rr);
%
% NOTE: All results are exploratory. Do not use for clinical decisions.

    if nargin < 2, r2_threshold = 0.93; end

    dfa_res = dfa(rr_ms);
    [chi, domain, interp] = chi_calculation(dfa_res.alpha1, dfa_res.alpha2);
    ph = phase_classification(dfa_res.alpha1, dfa_res.r_squared, chi, r2_threshold);

    result.alpha1             = dfa_res.alpha1;
    result.alpha2             = dfa_res.alpha2;
    result.r_squared          = dfa_res.r_squared;
    result.chi                = chi;
    result.chi_domain         = domain{1};
    result.chi_interpretation = interp{1};
    result.phase              = ph.label;
    result.phase_description  = ph.description;
    result.dfa                = dfa_res;

    fprintf('\n=== ECSoC Result [EXPLORATORY] ===\n');
    fprintf('  alpha1 = %.4f | alpha2 = %.4f | R2 = %.4f\n', result.alpha1, result.alpha2, result.r_squared);
    fprintf('  CHI    = %+.4f (%s)\n', result.chi, result.chi_domain);
    fprintf('  Phase  : %s\n', result.phase);
    fprintf('  %s\n', result.phase_description);
    fprintf('[Prospective validation required before any clinical use]\n\n');
end
