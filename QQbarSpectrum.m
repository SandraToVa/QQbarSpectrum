% Quarkonium spectrum with new potentials
% Different energies for states S, P, and D

%massa
%mass = load("dades.mat","m_c","m_b");
m_q = 1.496;

% arrays for the spectrum all can contain mixing

% ---------------------------------------------------------
% Calculation
% ---------------------------------------------------------

% Quarkonium spin 0 spectrum (orederd by J=L)
[s,ws,x0] = QQbarS0J0(m_q);
[p,wp,x1] = QQbarS0J1(m_q);
[d,wd,x2] = QQbarS0J2(m_q);

% Quarkonium spin 1 spectrum (ordered by J \neq L)
% Without mixing this gives the same as spin 0 states because we don't have
% hyperfine splitting of quarkonium at this order.
% We can see which state is it by the shape of the wave function
[j0,wj0,jx0] = QQbarS1J0(m_q);
[j1,wj1,jx1] = QQbarS1J1(m_q);
[j2,wj2,jx2] = QQbarS1J2(m_q);

% Hybrids spin 0 spectrum (ordered by Jcal=J \neq L)
% We can see which state is it by the shape of the wave function
[h0,wh0,hx0] = GQQbarS0Jcal0(m_q);
[h1,wh1,hx1] = GQQbarS0Jcal1(m_q);
[h2,wh2,hx2] = GQQbarS0Jcal2(m_q);

% Hybrids spin 1 spectrum (ordered by Jcal \neq J \neq L)
% Without mixing this gives the different as spin 0 states because we have
% hyperfine splitting of hybrids at this order.
% We can see which state is it by the shape of the wave function
A=-0.070;
B=0.0117;

[jh0,wjh0,jhx0] = GQQbarS1Jcal0(m_q,A,B);
[jh1,wjh1,jhx1] = GQQbarS1Jcal1(m_q,A,B);
[jh2,wjh2,jhx2] = GQQbarS1Jcal2(m_q,A,B);
[jh3,wjh3,jhx3] = GQQbarS1Jcal3(m_q,A,B);

% ---------------------------------------------------------
% Plots
% ---------------------------------------------------------

spin=true;
G=true;

if spin
    if G
        wf_all = {wjh0, wjh1, wjh2, wjh3};
        x_all  = {jhx0, jhx1, jhx2, jhx3};
        e_all  = {jh0, jh1, jh2, jh3};
        j_list = [0, 1, 2, 3];
    else
        wf_all = {wj0, wj1, wj2};
        x_all  = {jx0, jx1, jx2};
        e_all  = {j0, j1, j2};
        j_list = [0, 1, 2];
    end
else
    if G
        wf_all = {wh0, wh1, wh2};
        x_all  = {hx0, hx1, hx2};
        e_all  = {h0, h1, h2};
        j_list = [0, 1, 2];
    else
        wf_all = {ws, wp, wd};
        x_all  = {x0, x1, x2};
        e_all  = {s, p, d};
        j_list = [0, 1, 2];
    end
end

for n = 1:length(j_list)
    % Open a separate figure window for each J_cal
    figure('Color', 'w', 'Name', sprintf('J_cal = %d', j_list(n)));
    x_curr = x_all{n};
    
    for i = 1:9
        subplot(3, 3, i); % 3x3 layout fits up to 9 subplots
        
        % Extract all channel components for state i
        u_curr = squeeze(wf_all{n}(:, :, i));

        % Cycle through 7 colors, switching line style every 7 components
        % so every (color, style) pair is unique
        colors  = lines(7);
        lstyles = {'-', '--', ':', '-.'};
        hold on;
        for c = 1:size(u_curr, 1)
            plot(x_curr, u_curr(c, :), ...
                'Color', colors(mod(c-1, 7) + 1, :), ...
                'LineStyle', lstyles{mod(floor((c-1)/7), numel(lstyles)) + 1}, ...
                'LineWidth', 1.2);
        end
        hold off;
        grid on;
        xlim([0, 6]);
        xlabel('r [GeV^{-1}]');
        ylabel('u_c(r)');
        title(sprintf('State %d (E = %.4f GeV)', i, e_all{n}(i)));
  
    end
    % Dynamic Legend for channel components
    num_comp = size(u_curr, 1);
    comp_labels = arrayfun(@(c) sprintf('Comp %d', c), 1:num_comp, 'UniformOutput', false);
    legend(comp_labels, 'Location', 'best', 'FontSize', 7);
end

%%
% ---------------------------------------------------------
% Final hyperfine spectrum
% ---------------------------------------------------------
% Data
m_q = 1.496;
A_nomix=-0.0447;
B_nomix=0.0014;
A_mix=0;
B_mix=0;
% Order [{4 (s/d)1 states}, {4 p_1 states}, {4 (p/f)2 states}, {2 p_0 states}]
t = [4.0296 3.8976 3.9286 4.0746 4.1436 4.1106 4.1116 4.1756 ...
     4.2306 4.1786 4.2386 4.2516 4.4396 4.5136];

e = [0.0176 0.0186 0.0236 0.0216 0.0256 0.0276 0.0236 0.0186 ...
     0.0326 0.0276 0.0266 0.0346 0.0466 0.0536];

% 1. Compute spectrum
E_calc_nomix = compute_spectrum_nomix(m_q, A, B); 
E_calc_mix = compute_spectrum_mix(m_q, A, B);

% 2. Calculate chi-squared components
residuals_nomix = t - E_calc_nomix;             % Differences between exp (t) and calc
weighted_residuals_nomix = residuals_nomix ./ e; % Scaled by uncertainties (e)
chi2_nomix = sum(weighted_residuals_nomix.^2);   % Total chi^2
residuals_mix = t - E_calc_mix;             
weighted_residuals_mix = residuals_mix ./ e; 
chi2_mix = sum(weighted_residuals_mix.^2);   

% 3. Calculate degrees of freedom
N = length(t);   % Number of data points (e.g., 14)
p = 2;           % Number of fitted parameters (A and B)
dof = N - p;

% 4. Compute reduced chi-squared
chi2_red_nomix = chi2_nomix / dof;
chi2_red_mix = chi2_mix / dof;


% ---------------------------------------------------------
% Functions
% ---------------------------------------------------------

function E_list = compute_spectrum_nomix(m_q, A, B)
% Function to compute the spectrum without mixing

    % Hybrids spin 0 spectrum (ordered by Jcal=J \neq L)
    [h0,~,~] = GQQbarS0Jcal0(m_q);
    [h1,~,~] = GQQbarS0Jcal1(m_q);
    [h2,~,~] = GQQbarS0Jcal2(m_q);
    
    % Hybrids spin 1 spectrum
    % Compute the spectra with the current test values for A and B
    [jh0,~,~] = GQQbarS1Jcal0(m_q, A, B);
    [jh1,~,~] = GQQbarS1Jcal1(m_q, A, B);
    [jh2,~,~] = GQQbarS1Jcal2(m_q, A, B);
    [jh3,~,~] = GQQbarS1Jcal3(m_q, A, B);
    
    % Initialize the calculated energy vector
    a = zeros(1, 14);
    
    % MAP YOUR ENERGIES HERE
    % Assigning specific values of jhn to the vector 'a'
    a(1)  = h1(1); 
    a(2)  = jh0(1); 
    a(3)  = jh1(1); 
    a(4)  = jh2(1); 
    a(5)  = h1(2);
    a(6)  = jh0(2);
    a(7)  = jh1(2);
    a(8)  = jh2(2);
    a(9)  = h2(1);
    a(10) = jh1(3);
    a(11) = jh2(4);
    a(12) = jh3(1);
    a(13) = h0(1);
    a(14) = jh1(6);
    
    E_list = a;
end

function E_list = compute_spectrum_mix(m_q, A, B)
% Function to compute the spectrum with mixing

    % Hybrids spin 0 spectrum (ordered by Jcal=J \neq L)
    [h0,~,~] = GQQbarS0Jcal0(m_q);
    [h1,~,~] = GQQbarS0Jcal1(m_q);
    [h2,~,~] = GQQbarS0Jcal2(m_q);
    
    % Hybrids spin 1 spectrum
    % Compute the spectra with the current test values for A and B
    [jh0,~,~] = GQQbarS1Jcal0(m_q, A, B);
    [jh1,~,~] = GQQbarS1Jcal1(m_q, A, B);
    [jh2,~,~] = GQQbarS1Jcal2(m_q, A, B);
    [jh3,~,~] = GQQbarS1Jcal3(m_q, A, B);
    
    % Initialize the calculated energy vector
    a = zeros(1, 14);
    
    % MAP YOUR ENERGIES HERE
    % Assigning specific values of jhn to the vector 'a'
    a(1)  = h1(6); 
    a(2)  = jh0(3); 
    a(3)  = jh1(3); 
    a(4)  = jh2(3); 
    a(5)  = h1(9);
    a(6)  = jh0(4);
    a(7)  = jh1(4);
    a(8)  = jh2(2);
    a(9)  = h2(6);
    a(10) = jh1(5);
    a(11) = jh2(6);
    a(12) = jh3(2);
    a(13) = h0(4);
    a(14) = jh1(9);
    
    E_list = a;
end