% This code finds the best value for A and B (parameters of the
% interpolation of the hyperfine splitting) comparing with the charmonium
% hybrid lattice spectrum t (errors e).
%
% How it works
%  - The data t, e and the state assigned to each data point are read from
%    hybridAssignments.m (one table for mix=false and one for mix=true).
%  - The Hamiltonian of each spin-1 hybrid family (Jcal = 0,1,2,3) is exactly
%    linear in A and B:  H(A,B) = H0 + A*HA + B*HB. H0, HA, HB are built once
%    on the Lagrange-Laguerre mesh, so every (A,B) evaluation only needs a
%    diagonalization (no potential evaluations, no shooting, no NaNs).
%  - dE/dA = <psi|HA|psi>/m and dE/dB = <psi|HB|psi>/m (Hellmann-Feynman), so
%    lsqnonlin gets the exact Jacobian.
%  - The index in the assignment table is the position in energy order at
%    the best fit. The grid scan uses plain energy order. For the local fit
%    each level is followed as "k-th state of its decoupled channel block"
%    (so levels of different blocks can cross without being swapped, which
%    keeps chi^2 smooth). At the minimum the energy order is checked again
%    and, if it changed, the fit is repeated from there until it is stable.
%  - Errors from the Jacobian at the minimum: cov = inv(J'*J) (Delta chi^2=1),
%    also given scaled by sqrt(chi^2/dof) when chi^2/dof > 1 (PDG convention).
%    The chi^2 map is plotted with the 68% / 95% contours for 2 parameters.
%  - Errors of the fitted energies: the (A,B) covariance (with correlation)
%    is propagated to each E, checked on the 1-sigma ellipse, and added in
%    quadrature to sigma_th = 0.03 GeV (higher orders).

% =========================================================
% Settings
% =========================================================
clear; clc;
tic;

% Flags for all hybrid files (see hybridConfig.m)
hybridConfig('reset');
%cfg = hybridConfig('Eg', true, 'mix', false, 'hf', true);
cfg = hybridConfig('Eg', true, 'mix', true, 'hf', true);

% Quark mass
m_q = 1.496;

% Lattice data and assignment of each data point (see hybridAssignments.m)
[assign, t, e] = hybridAssignments(cfg.mix);

% Parameter bounds [A, B]
lb = [-0.3, -0.06];
ub = [ 0.3,  0.06];

% Solver mesh (N points per channel, last point at rmax GeV^-1).
% N=60 already gives the energies to ~1e-6 GeV (checked at the end)
solverOpts.N = 60;
solverOpts.rmax = 25;

% Grid used for the starting point and the chi^2 map
N_grid = 25;

% =========================================================
% 1. Build the model
% =========================================================
fprintf('Building Hamiltonians (mix=%d, hf=%d, Eg=%d)...\n', cfg.mix, cfg.hf, cfg.Eg);
model = build_model(m_q, assign, solverOpts);

% =========================================================
% 2. Coarse grid (starting point + chi^2 map), energy order
% =========================================================
A_vec = linspace(lb(1), ub(1), N_grid);
B_vec = linspace(lb(2), ub(2), N_grid);
chi2_map = zeros(N_grid, N_grid);   % chi2_map(j,i) <-> (A_vec(i), B_vec(j))
fprintf('Grid scan (%d x %d)...\n', N_grid, N_grid);
tgrid = tic;
for i = 1:N_grid
    for j = 1:N_grid
        E = model_spectrum(identify_states(model, A_vec(i), B_vec(j)), A_vec(i), B_vec(j));
        chi2_map(j,i) = sum(((t - E)./e).^2);
    end
end
fprintf('Grid scan done in %.1f s\n', toc(tgrid));
[chi2_grid, imin] = min(chi2_map(:));
[jmin, imin] = ind2sub(size(chi2_map), imin);
p0 = [A_vec(imin), B_vec(jmin)];
fprintf('Best grid point: A = %.4f, B = %.4f (chi^2 = %.3f)\n\n', p0(1), p0(2), chi2_grid);

% =========================================================
% 3. Local minimization (lsqnonlin, exact Jacobian), repeated until the
%    energy order at the minimum matches the one used in the fit
% =========================================================
options = optimoptions('lsqnonlin', ...
    'Display', 'final', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10);
max_refits = 5;
consistent = false;
for it = 1:max_refits
    model = identify_states(model, p0(1), p0(2));
    obj_fun = @(p) residuals(p, model, t, e);
    [best_params, chi2_min, ~, exitflag] = lsqnonlin(obj_fun, p0, lb, ub, options);
    model_best = identify_states(model, best_params(1), best_params(2));
    if isequal(model_best.block, model.block) && isequal(model_best.ord, model.ord)
        consistent = true;
        break;
    end
    fprintf('Energy order at the minimum differs from the start point: refitting (%d)\n', it);
    p0 = best_params;
end
if ~consistent
    warning(['The energy order at the minimum keeps changing: some assigned levels ' ...
             'are (nearly) degenerate with another level. Check the assignment table.']);
end
A_central = best_params(1);
B_central = best_params(2);

% =========================================================
% 4. Errors
% =========================================================
[~, J] = residuals(best_params, model, t, e);
cov_matrix = inv(J' * J);
param_errors = sqrt(diag(cov_matrix))';
rho_AB = cov_matrix(1,2) / prod(param_errors);

N = length(t);
p = 2;
dof = N - p;
chi2_red = chi2_min / dof;
scale = max(1, sqrt(chi2_red));

at_bound = abs(best_params - lb) < 1e-6*(ub-lb) | abs(best_params - ub) < 1e-6*(ub-lb);

% Mesh convergence check at the minimum
[E_central, dEdA, dEdB] = model_spectrum(model, A_central, B_central);
opts2 = solverOpts; opts2.N = round(1.4*solverOpts.N); opts2.rmax = 1.2*solverOpts.rmax;
model2 = identify_states(build_model(m_q, assign, opts2), A_central, B_central);
dE_mesh = max(abs(model_spectrum(model2, A_central, B_central) - E_central));

% =========================================================
% 5. Errors of the fitted energies
% =========================================================
% Covariance of (A,B) propagated to each energy, keeping the correlation:
%   sigma_fit^2 = G*C*G',  G = [dE/dA, dE/dB] (Hellmann-Feynman).
% Cross-check beyond linear order: largest |E - E_central| on the 1-sigma
% ellipse of (A,B) (for a linear E both give exactly the same number).
% The error used is the ellipse one; spin-0 states do not depend on A, B.
% Total error: sigma_fit added in quadrature to sigma_th (higher orders).
C_E = scale^2 * cov_matrix;          % scaled by chi^2/dof when > 1 (PDG)
sigma_th = 0.03;                     % GeV
sigma_lin = sqrt(sum(([dEdA(:), dEdB(:)] * C_E) .* [dEdA(:), dEdB(:)], 2))';
n_ell = 72;
Lc = chol(C_E, 'lower');
sigma_fit = zeros(1, N);
for q = 1:n_ell
    dp = Lc * [cos(2*pi*q/n_ell); sin(2*pi*q/n_ell)];
    Eq = model_spectrum(model, A_central + dp(1), B_central + dp(2));
    sigma_fit = max(sigma_fit, abs(Eq - E_central));
end
sigma_E = sqrt(sigma_fit.^2 + sigma_th^2);
if any(abs(sigma_fit - sigma_lin) > 0.1*sigma_lin + 1e-6)
    warning('Some energies are not linear in (A,B) over the 1-sigma ellipse: check sigma_lin vs sigma_fit.');
end

elapsed_time = toc;

print_assignment(model, t, A_central, B_central);

fprintf('\n================ Optimization Results ================\n');
fprintf('A = %.6f +/- %.6f   (scaled: +/- %.6f)\n', A_central, param_errors(1), scale*param_errors(1));
fprintf('B = %.6f +/- %.6f   (scaled: +/- %.6f)\n', B_central, param_errors(2), scale*param_errors(2));
fprintf('Correlation rho(A,B) = %.3f\n', rho_AB);
fprintf('------------------------------------------------------\n');
fprintf('Total Chi^2:        %.4f\n', chi2_min);
fprintf('Degrees of Freedom: %d\n', dof);
fprintf('Chi^2 / dof:        %.4f\n', chi2_red);
fprintf('p-value:            %.4f\n', 1 - chi2cdf(chi2_min, dof));
fprintf('lsqnonlin exitflag: %d, energy order consistent: %d (fits: %d)\n', exitflag, consistent, it);
fprintf('Mesh check: max |dE| (N=%d vs N=%d) = %.1e GeV\n', solverOpts.N, opts2.N, dE_mesh);
fprintf('Total Runtime:      %.2f seconds\n', elapsed_time);
fprintf('======================================================\n');
if any(at_bound)
    pnames = {'A','B'};
    warning('The minimum is at the boundary of [lb, ub] for %s: enlarge the bounds (errors are not reliable).', ...
        strjoin(pnames(at_bound), ' and '));
end
fprintf('\nErrors of E_fit: (A,B) covariance%s on the 1-sigma ellipse, total = sqrt(sigma_fit^2 + %.3f^2)\n', ...
    repmat(' (scaled)', 1, scale > 1), sigma_th);
fprintf('  #   data +/- err        E_fit +/- total    sigma_fit (lin)    pull  pull_tot\n');
for k = 1:N
    fprintf(' %2d  %.4f +/- %.4f   %.4f +/- %.4f   %.4f (%.4f)  %+6.2f  %+6.2f\n', k, t(k), e(k), ...
        E_central(k), sigma_E(k), sigma_fit(k), sigma_lin(k), (t(k)-E_central(k))/e(k), ...
        (t(k)-E_central(k))/sqrt(e(k)^2 + sigma_E(k)^2));
end

% =========================================================
% Plots
% =========================================================
figure('Color','w');
errorbar(1:N, t, e, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Lattice');
hold on;
errorbar((1:N) + 0.15, E_central, sigma_E, 'r.', 'MarkerSize', 18, 'CapSize', 3, ...
     'DisplayName', sprintf('Fit (\\chi^2/dof = %.2f), \\sigma = (\\sigma_{fit}^2 + %.2f^2)^{1/2}', chi2_red, sigma_th));
errorbar((1:N) + 0.15, E_central, sigma_fit, 'r', 'LineStyle', 'none', 'LineWidth', 2.5, 'CapSize', 0, ...
     'DisplayName', '\sigma_{fit} (A, B)');
xlabel('State index'); ylabel('Energy (GeV)');
title(sprintf('Charmonium hybrid spectrum fit: A = %.4f, B = %.4f', A_central, B_central));
legend('Location', 'best'); grid on;

figure('Color','w');
contourf(A_vec, B_vec, log10(chi2_map - min(chi2_map(:)) + 1), 30, 'LineColor', 'none');
colorbar; hold on;
% 68% and 95% regions for 2 parameters, from a fine grid around the minimum
A_fine = A_central + scale*param_errors(1)*linspace(-3.5, 3.5, 31);
B_fine = B_central + scale*param_errors(2)*linspace(-3.5, 3.5, 31);
chi2_fine = zeros(numel(B_fine), numel(A_fine));
for i = 1:numel(A_fine)
    for j = 1:numel(B_fine)
        chi2_fine(j,i) = sum(((t - model_spectrum(model, A_fine(i), B_fine(j)))./e).^2);
    end
end
contour(A_fine, B_fine, chi2_fine - chi2_min, [2.30 6.18], 'w', 'LineWidth', 1.5);
plot(A_central, B_central, 'w+', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('A'); ylabel('B');
title('log_{10}(\chi^2 - \chi^2_{min} + 1), contours \Delta\chi^2 = 2.30, 6.18');


% =========================================================
% Local Functions
% =========================================================

function model = build_model(m_q, assign, solverOpts)
% Precomputes everything that does not depend on A and B
    mesh = radialMesh(solverOpts.N, solverOpts.rmax);
    model.m = m_q;
    model.nData = size(assign, 1);
    model.spin = cell2mat(assign(:,1))';
    model.Jc = cell2mat(assign(:,2))';
    model.idx = cell2mat(assign(:,3))';

    % Spin-0 hybrids: independent of A, B
    E0 = cell(1,3);
    for Jc = 0:2
        f = str2func(sprintf('GQQbarS0Jcal%d', Jc));
        [~, ~, ~, info] = f(m_q, 'potential');
        H = radialHamiltonian(info.V, mesh);
        E0{Jc+1} = sort(eig(H))' / m_q;
    end

    % Spin-1 hybrids: H = H0 + A*HA + B*HB, split into decoupled blocks
    fam = cell(1,4);
    for Jc = 0:3
        f = str2func(sprintf('GQQbarS1Jcal%d', Jc));
        [~,~,~,i00] = f(m_q, 0, 0, 'potential');
        [~,~,~,i10] = f(m_q, 1, 0, 'potential');
        [~,~,~,i01] = f(m_q, 0, 1, 'potential');
        V00 = i00.V(mesh.r'); V10 = i10.V(mesh.r'); V01 = i01.V(mesh.r');
        n = size(V00, 1);
        VA = V10 - V00;  VB = V01 - V00;

        % Check the linearity in A and B at an arbitrary point
        Atest = -0.137; Btest = 0.0231;
        [~,~,~,itest] = f(m_q, Atest, Btest, 'potential');
        Vtest = itest.V(mesh.r');
        lin_err = max(abs(Vtest - (V00 + Atest*VA + Btest*VB)), [], 'all') / max(abs(VA), [], 'all');
        if lin_err > 1e-8
            error('The potential of GQQbarS1Jcal%d is not linear in A,B (err = %g).', Jc, lin_err);
        end

        % Decoupled channel blocks
        coupling = sum(abs(V00) + abs(VA) + abs(VB), 3) > 0;
        blk = conncomp(graph(coupling | coupling' | eye(n)));
        F.n = n;
        F.blocks = arrayfun(@(b) find(blk == b), 1:max(blk), 'UniformOutput', false);
        F.H0 = {}; F.HA = {}; F.HB = {};
        for b = 1:numel(F.blocks)
            ch = F.blocks{b};
            F.H0{b} = assemble(mesh, V00(ch,ch,:), true);
            F.HA{b} = assemble(mesh, VA(ch,ch,:), false);
            F.HB{b} = assemble(mesh, VB(ch,ch,:), false);
        end
        fam{Jc+1} = F;
    end

    model.E0 = E0;
    model.fam = fam;
    model.block = zeros(1, model.nData);
    model.ord = zeros(1, model.nData);
    model.mesh = mesh;
end


function model = identify_states(model, A, B)
% Finds, for each assigned spin-1 state (Jcal, index in energy order at
% (A,B)), its channel block and its position inside that block
    for Jc = 0:3
        ks = find(model.spin == 1 & model.Jc == Jc);
        if isempty(ks), continue; end
        F = model.fam{Jc+1};
        Eall = []; ball = []; oall = [];
        for b = 1:numel(F.blocks)
            Eb = sort(eig(F.H0{b} + A*F.HA{b} + B*F.HB{b}))';
            Eall = [Eall, Eb]; ball = [ball, b*ones(size(Eb))]; oall = [oall, 1:numel(Eb)]; %#ok<AGROW>
        end
        [~, s] = sort(Eall);
        for k = ks
            model.block(k) = ball(s(model.idx(k)));
            model.ord(k) = oall(s(model.idx(k)));
        end
    end
end


function H = assemble(mesh, V, withKinetic)
    n = size(V, 1);
    Vblk = cell(1, mesh.N);
    for i = 1:mesh.N
        Vblk{i} = V(:,:,i);
    end
    H = blkdiag(Vblk{:});
    if withKinetic
        H = H + kron(mesh.T, eye(n));
    end
    H = (H + H') / 2;
end


function [E, dEdA, dEdB] = model_spectrum(model, A, B)
% Energies (GeV) of the assigned states (followed by block and position in
% block, see identify_states) and their derivatives w.r.t. A, B
    E = zeros(1, model.nData);
    dEdA = zeros(1, model.nData);
    dEdB = zeros(1, model.nData);
    m = model.m;
    for k = find(model.spin == 0)
        E(k) = model.E0{model.Jc(k)+1}(model.idx(k));
    end
    for Jc = 0:3
        F = model.fam{Jc+1};
        ks = find(model.spin == 1 & model.Jc == Jc);
        for b = unique(model.block(ks))
            H = F.H0{b} + A*F.HA{b} + B*F.HB{b};
            kb = ks(model.block(ks) == b);
            nmax = max(model.ord(kb));
            if nargout > 1
                [C, lam] = eig(H, 'vector');
                [lam, s] = sort(lam); C = C(:, s(1:nmax));
            else
                lam = sort(eig(H));
            end
            for k = kb
                o = model.ord(k);
                E(k) = lam(o) / m;
                if nargout > 1
                    c = C(:, o);
                    dEdA(k) = (c' * F.HA{b} * c) / m;
                    dEdB(k) = (c' * F.HB{b} * c) / m;
                end
            end
        end
    end
end


function [res, J] = residuals(params, model, t, e)
    if nargout > 1
        [E, dA, dB] = model_spectrum(model, params(1), params(2));
        J = -[dA(:) ./ e(:), dB(:) ./ e(:)];
    else
        E = model_spectrum(model, params(1), params(2));
    end
    res = (t - E) ./ e;
end


function print_assignment(model, t, A, B)
% Table to check that each data point is matched to the intended state
    E = model_spectrum(model, A, B);
    fprintf('\nAssignment at the best fit A = %.4f, B = %.4f\n', A, B);
    fprintf('  #  spin Jcal idx  block(channels) #pos   E_fit    data\n');
    for k = 1:model.nData
        if model.spin(k) == 0
            bstr = '-';
        else
            ch = model.fam{model.Jc(k)+1}.blocks{model.block(k)};
            bstr = sprintf('%d(%s) #%d', model.block(k), strjoin(string(ch), ','), model.ord(k));
        end
        fprintf(' %2d   %d    %d   %d   %-20s %.4f  %.4f\n', k, model.spin(k), ...
            model.Jc(k), model.idx(k), bstr, E(k), t(k));
    end
    for Jc = 0:3
        if isscalar(model.fam{Jc+1}.blocks)
            fprintf('  Note: all channels of Jcal=%d are coupled, states are followed by energy order.\n', Jc);
        end
    end
end
