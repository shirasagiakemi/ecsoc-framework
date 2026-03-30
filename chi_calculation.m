function [chi, domain, interpretation] = chi_calculation(alpha1, alpha2)
% CHI_CALCULATION  Compute Criticality Heterogeneity Index (CHI)
%
% ECSoC Framework — MATLAB implementation
%
% SYNTAX:
%   [chi, domain] = chi_calculation(alpha1, alpha2)
%   [chi, domain, interpretation] = chi_calculation(alpha1, alpha2)
%
% FORMULA:
%   CHI = 2 × (α₁ − α₂)
%
% DOMAIN:
%   CHI >  0.05  → 'supercritical'   (dominant short-scale fluctuations)
%   CHI between -0.05 and 0.05 → 'near-critical'
%   CHI < -0.05  → 'subcritical'     (suppressed short-scale dynamics)
%
% REVISED REFERENCE FRAME (ECSoC v2, this paper):
%   Healthy (NSR, N=54):    mean CHI = +0.440, 87% supercritical
%   CHF NYHA1-3 (N=29):     mean CHI = −0.243, 28% supercritical
%   CAST baseline (N=734):  mean CHI = +0.113, 65.7% supercritical
%
% INPUTS:
%   alpha1 — Short-scale DFA exponent (scalar or vector)
%   alpha2 — Long-scale DFA exponent  (scalar or vector)
%
% OUTPUTS:
%   chi            — CHI value(s)
%   domain         — Cell array of domain strings
%   interpretation — Cell array of clinical interpretations [optional]

    chi = 2 .* (alpha1 - alpha2);

    % Domain classification
    domain = cell(size(chi));
    for i = 1:numel(chi)
        if chi(i) > 0.05
            domain{i} = 'supercritical';
        elseif chi(i) < -0.05
            domain{i} = 'subcritical';
        else
            domain{i} = 'near-critical';
        end
    end

    if nargout < 3, return; end

    % Interpretation relative to healthy reference
    healthy_mean = 0.440;
    chf_mean     = -0.243;

    interpretation = cell(size(chi));
    for i = 1:numel(chi)
        if chi(i) >= healthy_mean - 0.441
            interpretation{i} = 'Within healthy supercritical range';
        elseif chi(i) >= 0
            interpretation{i} = 'Supercritical but below healthy norm';
        elseif chi(i) >= chf_mean
            interpretation{i} = 'Subcritical shift; consistent with mild-moderate CHF';
        else
            interpretation{i} = 'Markedly subcritical; consistent with advanced disease';
        end
    end
end
