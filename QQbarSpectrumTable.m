% Spectrum tables in the format of the tables of mixing.pdf (Figs. 5 and 9)
% For every state: matrix, n, S, J_tot, j, L, E, J^PC, composition
% (dominant channel), quarkonium/hybrid and %Q, %H.
%
% n    : position of the state in the energy order of its file
% S    : quark spin, j : J of the dominant channel (gluon + L for hybrids,
%        J_tot for quarkonium), L : orbital angular momentum
% E    : Schroedinger energy (MeV), without the 2m - Eval shift
% M    : mass of the state, E + 2m - Eval (MeV)
% P = (-1)^(L+1), C = cg (-1)^(L+S), cg = +1 quarkonium, -1 hybrid,
% with L and S of the dominant channel.
%
% Autor: Sandra Tomàs

m_q = 1.496;
A = -0.044691;
B = 0.001377;

Emax = 1.80;   % highest E (GeV) listed
nmax = 20;     % states computed per file

hybridConfig('reset');
quarkoniumConfig('reset');
hybridConfig('mix', true, 'hf', true);
quarkoniumConfig('mix', true);

% Quarkonium files (as in Fig. 5)
quarkoniumFiles = {'QQbarS0J0','QQbarS0J1','QQbarS0J2', ...
                   'QQbarS1J0','QQbarS1J1','QQbarS1J2'};
% Hybrid files (as in Fig. 9)
hybridFiles = {'GQQbarS0Jcal0','GQQbarS0Jcal1','GQQbarS0Jcal2', ...
               'GQQbarS1Jcal0','GQQbarS1Jcal1','GQQbarS1Jcal2','GQQbarS1Jcal3'};

Tq = spectrumTable(quarkoniumFiles, m_q, A, B, Emax, nmax);
printSpectrumTable(Tq, 'Quarkonium files');

Th = spectrumTable(hybridFiles, m_q, A, B, Emax, nmax);
printSpectrumTable(Th, sprintf('Hybrid files (A = %g, B = %g)', A, B));

% writetable(Tq, 'spectrumQuarkonium.csv');
% writetable(Th, 'spectrumHybrid.csv');


% ---------------------------------------------------------
% Functions
% ---------------------------------------------------------

function T = spectrumTable(files, m, A, B, Emax, nmax)
% Solves every file and labels its states; T is sorted by energy
rows = {};
for k = 1:numel(files)
    name = files{k};
    J = str2double(name(end));
    isHybridFile = startsWith(name, 'G');
    if isHybridFile, cfg = hybridConfig; else, cfg = quarkoniumConfig; end

    % Potential of the file and its channels
    f = str2func(name);
    try
        if startsWith(name, 'GQQbarS1')
            [~,~,~,info] = f(m, A, B, 'potential');
        else
            [~,~,~,info] = f(m, 'potential');
        end
        V = info.V(1);
    catch err
        warning('%s skipped: %s', name, err.message);
        continue;
    end
    [ch, block] = fileChannels(name, J);
    if ~cfg.mix
        % without mixing the file keeps only its own kind of channels
        ch = ch([ch.isQ] ~= isHybridFile);
    end
    if numel(ch) ~= size(V,1)
        error('%s: %d channels in the potential, %d labels.', name, size(V,1), numel(ch));
    end

    [E, ~, ~, sol] = solveRadialSchrodinger(info.V, m, nmax, struct());
    shift = 2*m - cfg.Eval*(m == 1.496);   % the files add it only if Eg = true
    Es = E - cfg.Eg*shift;
    if Es(end) < Emax
        warning('%s: all %d states are below Emax, increase nmax.', name, nmax);
    end

    for s = find(Es <= Emax)
        w = sol.weights(:, s);
        [~, d] = max(w);
        c = ch(d);
        pQ = 100*sum(w([ch.isQ]));
        if c.isQ, type = 'Quarkonium'; cg = 1; else, type = 'Hybrid'; cg = -1; end
        P = (-1)^(c.L+1);
        C = cg*(-1)^(c.L+c.S);
        JPC = sprintf('%d^{%s%s}', J, pm(P), pm(C));
        rows(end+1,:) = {block, s, c.S, J, c.j, c.L, 1000*Es(s), 1000*(Es(s)+shift), ...
                         JPC, c.label, type, pQ, 100-pQ, name}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Matrix','n','S','J','j','L', ...
    'E','M','JPC','Composition','Type','pQ','pH','File'});
T = sortrows(T, 'E');
end


function [ch, block] = fileChannels(name, J)
% Channels of each file, in the order of its potential matrix (mix = true)
if startsWith(name, 'QQbarS0J') || startsWith(name, 'GQQbarS1Jcal')
    % MixingH1: R0 and the spin-1 hybrid P1^{ab}
    block = 'MixingH1';
    if J == 0
        ch = [R0(J), P1(J,1,1), P1(J,-1,1)];
    elseif J == 1
        ch = [R0(J), P1(J,1,1), P1(J,-1,1), P1(J,1,-1), P1(J,0,0)];
    else
        ch = [R0(J), P1(J,1,1), P1(J,-1,1), P1(J,1,-1), P1(J,-1,-1), P1(J,0,0)];
    end
    if startsWith(name, 'GQQbarS1Jcal')
        % spin-1 hybrid channels decoupled from quarkonium
        if J == 0
            ch = [ch, P1(J,0,1)];
        elseif J == 1
            ch = [ch, P1(J,0,1), P1(J,1,0), P1(J,-1,0)];
        else
            ch = [ch, P1(J,0,1), P1(J,0,-1), P1(J,1,0), P1(J,-1,0)];
        end
    end
elseif startsWith(name, 'QQbarS1J') || startsWith(name, 'GQQbarS0Jcal')
    % MixingH0: R1^a and the spin-0 hybrid P0^a
    block = 'MixingH0';
    if J == 0
        ch = [R1(J,1), P0(J,1)];
    else
        ch = [R1(J,0), P0(J,0), R1(J,1), R1(J,-1), P0(J,1), P0(J,-1)];
    end
else
    error('Unknown file %s.', name);
end
end

% a, b = -1, 0, +1 : L = J + a, J = Jcal + b (as in P^{ab} of the notes)
function c = R0(J),      c = channel('R0', 1, 0, J, J);              end
function c = R1(J,a),    c = channel(['R1^' pm(a)], 1, 1, J, J+a);   end
function c = P0(J,a),    c = channel(['P0^' pm(a)], 0, 0, J, J+a);   end
function c = P1(J,a,b),  c = channel(['P1^{' pm(a) pm(b) '}'], 0, 1, J+b, J+b+a); end

function c = channel(label, isQ, S, j, L)
c = struct('label', label, 'isQ', logical(isQ), 'S', S, 'j', j, 'L', L);
end

function s = pm(x)
s = '0';
if x > 0, s = '+'; elseif x < 0, s = '-'; end
end


function printSpectrumTable(T, titleText)
fprintf('\n%s\n', titleText);
fprintf('%-9s %3s %2s %4s %2s %2s %8s %8s %-7s %-9s %-11s %4s %4s  %s\n', ...
    'Matrix','n','S','JTot','j','L','E(MeV)','M(MeV)','JPC','Comp.','Type','%Q','%H','File');
for k = 1:height(T)
    fprintf('%-9s %3d %2d %4d %2d %2d %8.1f %8.1f %-7s %-9s %-11s %4.0f %4.0f  %s\n', ...
        T.Matrix{k}, T.n(k), T.S(k), T.J(k), T.j(k), T.L(k), T.E(k), T.M(k), ...
        T.JPC{k}, T.Composition{k}, T.Type{k}, T.pQ(k), T.pH(k), T.File{k});
end
end
