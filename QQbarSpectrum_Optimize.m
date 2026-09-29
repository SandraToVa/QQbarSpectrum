
% This code find the best value for A and B (parameters of the
% interpolation of the hyperfine splitting) comparing with charmonium
% experimental lattice data.

% =========================================================
% Hyperfine Splitting Fit (A, B Parameters) 
% =========================================================
clear; clc;
tic;

% Quark mass
m_q = 1.496;

% Experimental data
t = [4.0296 3.8976 3.9286 4.0746 4.1436 4.1106 4.1116 4.1756 ...
     4.2306 4.1786 4.2386 4.2516 4.4396 4.5136];
e = [0.0176 0.0186 0.0236 0.0216 0.0256 0.0276 0.0236 0.0186 ...
     0.0326 0.0276 0.0266 0.0346 0.0466 0.0536];

% Parameter Bounds [A, B]
lb = [-0.1, -0.06];  
ub = [ 0.1,  0.06];  

% 1. Precompute spin-0 spectrum
fprintf('Precomputing spin-0 states...\n');
h_fixed = precompute_spin0(m_q);

% 2. Coarse Grid Search
% 30x30=900 grid
N_grid = 30;
A_vec = linspace(lb(1), ub(1), N_grid);
B_vec = linspace(lb(2), ub(2), N_grid);
best_chi2 = Inf;
p0 = [0, 0];
skipped_points = 0;

fprintf('Starting Grid Search (%d x %d evaluations)...\n', N_grid, N_grid);
for i = 1:N_grid
    if mod(i, floor(N_grid/10)) == 0 || i == N_grid
        fprintf('  Grid search progress: %3d%%\n', round((i / N_grid) * 100));
    end
    
    for j = 1:N_grid
        [res, failed] = compute_residuals([A_vec(i), B_vec(j)], m_q, t, e, h_fixed);
        if failed
            skipped_points = skipped_points + 1;
            continue; % Skip invalid point
        end
        
        chi2_val = sum(res.^2);
        if chi2_val < best_chi2
            best_chi2 = chi2_val;
            p0 = [A_vec(i), B_vec(j)];
        end
    end
end

fprintf('Grid search finished. Skipped %d problematic points.\n', skipped_points);
fprintf('Best starting guess found: A = %.4f, B = %.4f (Initial Chi^2 = %.2f)\n\n', p0(1), p0(2), best_chi2);

% 3. Optimization using lsqnonlin with penalty boundary handling
fprintf('Starting lsqnonlin local optimization...\n');
options = optimoptions('lsqnonlin', ...
                       'Display', 'iter', ...
                       'FunctionTolerance', 1e-6, ...
                       'StepTolerance', 1e-6);

obj_fun = @(params) compute_residuals(params, m_q, t, e, h_fixed);
[best_params, resnorm, ~, ~, ~, ~, jacobian] = ...
    lsqnonlin(obj_fun, p0, lb, ub, options);

A_central = best_params(1);
B_central = best_params(2);

% 4. Statistical Error Estimation from Jacobian
J = full(jacobian);
cov_matrix = inv(J' * J); 
param_errors = sqrt(diag(cov_matrix));

err_A = param_errors(1);
err_B = param_errors(2);

% 5. Fit Quality Statistics
N = length(t);   
p = 2;            
dof = N - p;      
chi2_red = resnorm / dof;

elapsed_time = toc; 

% Display Results
fprintf('\n================ Optimization Results ================\n');
fprintf('A = %.6f +/- %.6f\n', A_central, err_A);
fprintf('B = %.6f +/- %.6f\n', B_central, err_B);
fprintf('------------------------------------------------------\n');
fprintf('Total Chi^2:        %.4f\n', resnorm);
fprintf('Degrees of Freedom: %d\n', dof);
fprintf('Chi^2 / dof:        %.4f\n', chi2_red);
fprintf('Total Runtime:      %.2f seconds\n', elapsed_time);
fprintf('======================================================\n');

% ---------------------------------------------------------
% Plot
% ---------------------------------------------------------
[E_calc, failed] = compute_spectrum(m_q, A_central, B_central, h_fixed);

if ~failed
    figure;
    errorbar(1:N, t, e, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Experimental');
    hold on;
    plot(1:N, E_calc, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12, ...
         'DisplayName', sprintf('Central Fit (\\chi^2/dof = %.2f)', chi2_red));
    xlabel('State index'); ylabel('Energy (GeV)');
    title('Charmonium Spectrum Fit');
    legend('Location', 'best');
    grid on;
else
    warning('Best parameters returned an invalid spectrum evaluation.');
end

% =========================================================
% Local Functions
% =========================================================

function h_fixed = precompute_spin0(m_q)
% Precompute the spectrum states where no A and B are involved
    [h0,~,~] = GQQbarS0Jcal0(m_q);
    [h1,~,~] = GQQbarS0Jcal1(m_q);
    [h2,~,~] = GQQbarS0Jcal2(m_q);
    h_fixed.h0 = h0;
    h_fixed.h1 = h1;
    h_fixed.h2 = h2;
end

function [E_calc, failed] = compute_spectrum(m_q, A, B, h_fixed)
    failed = false;
    E_calc = zeros(1, 14);

    % 1. Clear global MATLAB warning buffer before calculation
    lastwarn(''); 

    try
        % Execute physics calculation (calls GQQbar -> shotForEigenvalues)
        [jh0,~,~] = GQQbarS1Jcal0(m_q, A, B);
        [jh1,~,~] = GQQbarS1Jcal1(m_q, A, B);
        [jh2,~,~] = GQQbarS1Jcal2(m_q, A, B);
        [jh3,~,~] = GQQbarS1Jcal3(m_q, A, B);

        % 2. Check if shotForEigenvalues (or any other subfunction) issued a warning
        [last_msg, ~] = lastwarn;
        if ~isempty(last_msg)
            failed = true;
            return;
        end

        % 3. Check for non-physical / non-converged values
        if any(isnan(jh0)) || any(isnan(jh1)) || any(isnan(jh2)) || any(isnan(jh3)) || ...
           any(isinf(jh0)) || any(isinf(jh1)) || any(isinf(jh2)) || any(isinf(jh3))
            failed = true;
            return;
        end

        % Map spectrum array
        E_calc(1)  = h_fixed.h1(1); 
        E_calc(2)  = jh0(1); 
        E_calc(3)  = jh1(1); 
        E_calc(4)  = jh2(1); 
        E_calc(5)  = h_fixed.h1(2);
        E_calc(6)  = jh0(2);
        E_calc(7)  = jh1(2);
        E_calc(8)  = jh2(2);
        E_calc(9)  = h_fixed.h2(1);
        E_calc(10) = jh1(3);
        E_calc(11) = jh2(4);
        E_calc(12) = jh3(1);
        E_calc(13) = h_fixed.h0(1);
        E_calc(14) = jh1(6);

    catch
        % Triggers if a hard error occurs
        failed = true;
    end
end

function [res, failed] = compute_residuals(params, m_q, t, e, h_fixed)
    A = params(1);
    B = params(2);
    
    [E_calc, failed] = compute_spectrum(m_q, A, B, h_fixed);
    
    if failed
        % Assign a large penalty vector if calculation fails or emits warnings
        res = ones(size(t)) * 1e5;
    else
        res = (t - E_calc) ./ e;
    end
end