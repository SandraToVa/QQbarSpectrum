% Proposes the mix=true assignment table of hybridAssignments.m from the
% mix=false one.
%
% For every data point, the assigned mix=false state is compared with all
% the mix=true states of the same file (same mesh): the mix=false channels
% are located inside the mix=true matrix (their diagonal potentials are the
% same), and the mix=true state with the largest overlap
% |<psi_mix=false|psi_mix=true>| is proposed. The table also prints the
% overlap, the second best overlap and the quarkonium fraction (weight in
% the channels that only exist with mixing), so ambiguous cases are visible.
% Copy the printed table into hybridAssignments.m (mix=true case).

clear; clc;

% Point where the states are compared (e.g. the mix=false best fit)
m_q = 1.496;
A = -0.0548;
B =  0.0038;

nstates = 25;                 % states computed per file with mixing
opts.N = 60; opts.rmax = 25;  % same mesh for both cases

hybridConfig('reset');
cfg0 = hybridConfig;
assignF = hybridAssignments(false);
mesh = radialMesh(opts.N, opts.rmax);

nData = size(assignF, 1);
proposal = assignF;
fprintf('  #  spin Jcal | idx(mix=false) -> idx(mix=true)  overlap  2nd best  quarkonium frac   E_false  E_true\n');
for k = 1:nData
    spin = assignF{k,1}; Jc = assignF{k,2}; idxF = assignF{k,3};
    if spin == 0
        f = str2func(sprintf('GQQbarS0Jcal%d', Jc)); args = {m_q};
    else
        f = str2func(sprintf('GQQbarS1Jcal%d', Jc)); args = {m_q, A, B};
    end

    % The potential reads the mix flag when it is evaluated, so evaluate it
    % right after setting the flag
    hybridConfig('mix', false);
    [~,~,~,iF] = f(args{:}, 'potential');
    VF = iF.V(mesh.r');
    hybridConfig('mix', true);
    [~,~,~,iT] = f(args{:}, 'potential');
    VT = iT.V(mesh.r');
    hybridConfig('mix', cfg0.mix);

    [EF,~,~,infoF] = solveRadialSchrodinger([], m_q, idxF, struct('H', radialHamiltonian(@(r) VF, mesh), 'mesh', mesh));
    [ET,~,~,infoT] = solveRadialSchrodinger([], m_q, nstates, struct('H', radialHamiltonian(@(r) VT, mesh), 'mesh', mesh));
    nF = size(VF,1); nT = size(VT,1);

    % Locate each mix=false channel in the mix=true matrix
    map = zeros(1, nF);
    for c = 1:nF
        d = squeeze(VF(c,c,:));
        for c2 = 1:nT
            if max(abs(d - squeeze(VT(c2,c2,:)))) <= 1e-10*max(abs(d))
                map(c) = c2; break;
            end
        end
    end
    if any(map == 0)
        % The mix=false channels are not all present with mixing: no overlap
        % can be computed. Print the channel content of the mix=true states.
        proposal{k,3} = NaN;
        fprintf(' %2d   %d    %d  |      %d         ->      ??   channel %s of mix=false not found with mixing\n', ...
            k, spin, Jc, idxF, strjoin(string(find(map == 0)), ','));
        fprintf('       mix=true states: E and channel weights\n');
        for s = 1:min(6, nstates)
            fprintf('         %2d  %.4f   %s\n', s, ET(s), sprintf('%.3f ', infoT.weights(:,s)));
        end
        continue;
    end
    extra = setdiff(1:nT, map);   % channels that only exist with mixing

    % Overlaps with all mix=true states
    cF = reshape(infoF.C(:, idxF), nF, opts.N);
    cFT = zeros(nT, opts.N);
    cFT(map,:) = cF;
    ov = abs(cFT(:)' * infoT.C);
    [ovs, order] = sort(ov, 'descend');
    idxT = order(1);
    qfrac = sum(infoT.weights(extra, idxT));
    proposal{k,3} = idxT;
    fprintf(' %2d   %d    %d  |      %d         ->      %2d          %.3f    %.3f        %.3f          %.4f  %.4f\n', ...
        k, spin, Jc, idxF, idxT, ovs(1), ovs(2), qfrac, EF(idxF), ET(idxT));
end

fprintf('\nProposed mix=true table (copy into hybridAssignments.m):\n');
fprintf('    assign = { ');
for k = 1:nData
    fprintf('%d,%d,%d', proposal{k,:});
    if k < nData
        if mod(k,4) == 0, fprintf('; ...\n               '); else, fprintf(';   '); end
    end
end
fprintf(' };\n');
