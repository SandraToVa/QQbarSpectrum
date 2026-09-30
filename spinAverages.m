% Spin averages of the charmonium hybrid multiplets compared with the
% lattice data, at given A and B.
%
% For each multiplet ((s/d)1, p1, (p/f)2, p0) the spin average is
%     M_avg = sum_J (2J+1) M_J / sum_J (2J+1)
% over its members (the spin-0 state and the spin-1 states). The
% spin-dependent (hyperfine) terms mostly cancel in this average, so a
% difference between model and data in M_avg points to the static
% potentials, while the splittings inside a multiplet test the hyperfine
% terms (A, B).
%
% The data errors are treated as uncorrelated, as in QQbarSpectrum_Optimize.
% The data, the assignment of each point and its J^PC are read from
% hybridAssignments.m.

clear; clc;

% =========================================================
% Settings
% =========================================================
hybridConfig('reset');
cfg_h = hybridConfig('Eg', true, 'mix', false, 'hf', true);

% Data
m_q = 1.496;
% A and B: best fit stored as ref in hybridAssignments.m
[assign, t, e, states, ref] = hybridAssignments(cfg_h.mix);
A = ref(1);
B = ref(2);

if cfg_h.mix
    disp('The following results are computing WITH mixing');
    disp('------------------------------------------------------');
else
    disp('The following results are computing WITHOUT mixing')
    disp('------------------------------------------------------');
end

nParams = 2;       % fitted parameters, for chi^2/dof
solverOpts.N = 60;
solverOpts.rmax = 25;

% =========================================================
% Model spectrum of the assigned states
% =========================================================
E = assigned_spectrum(assign, m_q, A, B, solverOpts);

nData = numel(t);
pull = (t - E) ./ e;
chi2 = sum(pull.^2);
dof = nData - nParams;
w = 2*states.J + 1;
names = states.multipletNames;
nMult = numel(names);

fprintf('Spin averages at A = %.6f, B = %.6f  (mix=%d, hf=%d, Eg=%d)\n\n', ...
    A, B, cfg_h.mix, cfg_h.hf, cfg_h.Eg);

% ---------------------------------------------------------
% Table 1: every data point
% ---------------------------------------------------------
fprintf('  #  multiplet  J^PC   data     model    diff(MeV)   pull\n');
for k = 1:nData
    fprintf(' %2d  %-9s  %-4s  %.4f   %.4f   %+7.1f    %+6.2f\n', k, ...
        names{states.multiplet(k)}, states.JPC{k}, t(k), E(k), 1000*(E(k)-t(k)), pull(k));
end
fprintf('\n');

% ---------------------------------------------------------
% Table 2: spin averages
% ---------------------------------------------------------
avg_data = zeros(1,nMult); avg_model = zeros(1,nMult); avg_err = zeros(1,nMult);
chi2_mult = zeros(1,nMult);
fprintf('multiplet   spin avg data   model    diff(MeV)   err(MeV)   pull  | chi2 of its points\n');
for q = 1:nMult
    ix = find(states.multiplet == q);
    W = w(ix) / sum(w(ix));
    avg_data(q)  = sum(W .* t(ix));
    avg_model(q) = sum(W .* E(ix));
    avg_err(q)   = sqrt(sum(W.^2 .* e(ix).^2));
    chi2_mult(q) = sum(pull(ix).^2);
    fprintf('%-9s   %.4f          %.4f   %+7.1f     %5.1f    %+5.2f  | %5.1f\n', names{q}, ...
        avg_data(q), avg_model(q), 1000*(avg_model(q)-avg_data(q)), 1000*avg_err(q), ...
        (avg_data(q)-avg_model(q))/avg_err(q), chi2_mult(q));
end
fprintf('\n');

% ---------------------------------------------------------
% Table 3: splittings inside each multiplet (w.r.t. its spin average)
% ---------------------------------------------------------
fprintf('Splittings M_J - M_avg (MeV), data / model\n');
for q = 1:nMult
    ix = find(states.multiplet == q);
    fprintf('  %-9s', names{q});
    for k = ix
        fprintf('  %s: %+5.0f / %+5.0f', states.JPC{k}, ...
            1000*(t(k)-avg_data(q)), 1000*(E(k)-avg_model(q)));
    end
    fprintf('\n');
end
fprintf('\n');

% ---------------------------------------------------------
% chi^2 breakdown
% ---------------------------------------------------------
% chi^2 of the spin averages alone
chi2_avg = sum(((avg_data - avg_model) ./ avg_err).^2);

% chi^2 left if each multiplet could be shifted freely (at fixed A, B)
chi2_shift = 0;
for q = 1:nMult
    ix = find(states.multiplet == q);
    d = t(ix) - E(ix); wt = 1 ./ e(ix).^2;
    s = sum(d .* wt) / sum(wt);
    chi2_shift = chi2_shift + sum(((d - s) ./ e(ix)).^2);
end

% best common shift of all the levels (e.g. a change of Eg)
d = t - E; wt = 1 ./ e.^2;
s_common = sum(d .* wt) / sum(wt);
chi2_common = sum(((d - s_common) ./ e).^2);

fprintf('==================== chi^2 breakdown ====================\n');
fprintf('Total chi^2                         %6.2f   chi^2/dof = %.2f  (dof %d)\n', chi2, chi2/dof, dof);
for q = 1:nMult
    fprintf('  from %-9s points              %6.2f\n', names{q}, chi2_mult(q));
end
fprintf('chi^2 of the %d spin averages alone   %6.2f\n', nMult, chi2_avg);
fprintf('Left if each multiplet shifts freely %6.2f   chi^2/dof = %.2f  (dof %d)\n', ...
    chi2_shift, chi2_shift/(dof-nMult), dof-nMult);
fprintf('Best common shift %+6.1f MeV         %6.2f   chi^2/dof = %.2f  (dof %d)\n', ...
    1000*s_common, chi2_common, chi2_common/(dof-1), dof-1);
fprintf('   (equivalent to Eval = %.4f GeV instead of %.4f)\n', cfg_h.Eval - s_common, cfg_h.Eval);
fprintf('Shifts are evaluated at fixed A, B (no refit).\n');
fprintf('=========================================================\n');


% =========================================================
% Local functions
% =========================================================

function E = assigned_spectrum(assign, m_q, A, B, solverOpts)
% Energies of the assigned states: index = position in energy order at (A,B)
    mesh = radialMesh(solverOpts.N, solverOpts.rmax);
    nData = size(assign, 1);
    E = zeros(1, nData);
    cache = containers.Map();
    for k = 1:nData
        spin = assign{k,1}; Jc = assign{k,2}; idx = assign{k,3};
        key = sprintf('%d_%d', spin, Jc);
        if ~isKey(cache, key)
            if spin == 0
                f = str2func(sprintf('GQQbarS0Jcal%d', Jc)); args = {m_q};
            else
                f = str2func(sprintf('GQQbarS1Jcal%d', Jc)); args = {m_q, A, B};
            end
            [~, ~, ~, info] = f(args{:}, 'potential');
            H = radialHamiltonian(info.V, mesh);
            cache(key) = sort(eig(H))' / m_q;
        end
        Ek = cache(key);
        E(k) = Ek(idx);
    end
end
