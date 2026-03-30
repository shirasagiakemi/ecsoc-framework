function phase = phase_classification(alpha1, r_squared, chi, r2_threshold)
% PHASE_CLASSIFICATION  ECSoC phase assignment (Phases I–V)
%
% ECSoC Framework — MATLAB implementation
%
% SYNTAX:
%   phase = phase_classification(alpha1, r_squared, chi)
%   phase = phase_classification(alpha1, r_squared, chi, r2_threshold)
%
% PHASE DEFINITIONS (thresholds from MUSIC cohort — exploratory):
%   Phase I   : α₁ < 0.85,  R² ≥ threshold  → Subcritical stable
%   Phase II  : 0.85 ≤ α₁ < 1.15, R² ≥ threshold → Near-critical stable
%   Phase III : α₁ ≥ 1.15, R² ≥ threshold   → Supercritical stable
%   Phase Va  : R² < threshold, CHI < -0.05  → Structural breakdown, subcritical
%   Phase Vb  : R² < threshold, CHI >  0.05  → Structural breakdown, supercritical
%   Phase Vc  : R² < threshold, |CHI| ≤ 0.05 → Structural breakdown, near-critical
%
% INPUTS:
%   alpha1       — Short-scale DFA exponent
%   r_squared    — DFA goodness-of-fit R²
%   chi          — CHI value (required for Phase V subtyping)
%   r2_threshold — R² cutoff for Phase V (default: 0.93)
%
% OUTPUT (struct):
%   .phase       — 'I', 'II', 'III', or 'V'
%   .subtype     — 'Va', 'Vb', 'Vc', or '' (for phases I–III)
%   .label       — Full label string (e.g. 'Phase Vb')
%   .description — Mechanistic interpretation string

    if nargin < 4, r2_threshold = 0.93; end

    if r_squared < r2_threshold
        % Phase V — structural breakdown
        if chi < -0.05
            phase.subtype     = 'Va';
            phase.description = 'Structural breakdown, subcritical: rigid low-variability with scaling loss';
        elseif chi > 0.05
            phase.subtype     = 'Vb';
            phase.description = 'Structural breakdown, supercritical: Path 2 collapse candidate';
        else
            phase.subtype     = 'Vc';
            phase.description = 'Structural breakdown, near-critical manifold';
        end
        phase.phase = 'V';
        phase.label = ['Phase ' phase.subtype];
        return
    end

    % Phases I–III
    phase.subtype = '';
    if alpha1 < 0.85
        phase.phase       = 'I';
        phase.description = 'Subcritical stable: rigid, suppressed autonomic modulation';
    elseif alpha1 < 1.15
        phase.phase       = 'II';
        phase.description = 'Near-critical stable: proximity to critical manifold (α₁ ≈ 1.0)';
    else
        phase.phase       = 'III';
        phase.description = 'Supercritical stable: dominant short-scale fluctuations (healthy norm)';
    end
    phase.label = ['Phase ' phase.phase];
end
