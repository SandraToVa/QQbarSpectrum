function [assign, t, e, states, ref] = hybridAssignments(mix)
% [assign, t, e, states, ref] = hybridAssignments(mix)
% Lattice data used in QQbarSpectrum_Optimize and the computed state that
% each data point is compared with. Row k of assign corresponds to t(k), e(k).
%
% assign{k,:} = {spin, Jcal, index}
%   spin  : 0 -> GQQbarS0Jcal<Jcal>.m,  1 -> GQQbarS1Jcal<Jcal>.m
%   index : position of the state in energy order (1 = lowest) within that
%           file, at the best-fit (A,B) = ref.
%
% ref = [A, B] where the indices are valid (the best fit). The optimizer
% identifies each spin-1 state at ref as "k-th level of its decoupled channel
% block" and follows that label during the whole fit (levels of the same
% block do not cross). If the energy order at the new minimum differs from
% the table, it prints the new indices: update the table and ref.
%
% If the potentials change: follow the steps "IF THE POTENTIALS CHANGE" at
% the top of QQbarSpectrum_Optimize.m (check the states at ref with the
% plots of QQbarSpectrum.m, refit, update indices + ref here and the same
% indices in compute_spectrum_nomix / compute_spectrum_mix of
% QQbarSpectrum.m).
%
% The indices change with mix (quarkonium states appear in the spectrum),
% so there is one table per case. proposeMixAssignment.m suggests the
% mix=true table from the mix=false one (overlap of wave functions).
%
% states (used by spinAverages.m) describes each data point:
%   states.JPC{k}      : J^PC label
%   states.J(k)        : total angular momentum J
%   states.multiplet(k): index into states.multipletNames

% Lattice data (order: {4 (s/d)1 states}, {4 p_1 states}, {4 (p/f)2 states},
% {2 p_0 states})
t = [4.0296 3.8976 3.9286 4.0746 4.1436 4.1106 4.1116 4.1756 ...
     4.2306 4.1786 4.2386 4.2516 4.4396 4.5136];
e = [0.0176 0.0186 0.0236 0.0216 0.0256 0.0276 0.0236 0.0186 ...
     0.0326 0.0276 0.0266 0.0346 0.0466 0.0536];

% Quantum numbers of each data point (same order as t)
states.JPC = {'1--','0-+','1-+','2-+', '1++','0+-','1+-','2+-', ...
              '2++','1+-','2+-','3+-', '0++','1+-'};
states.J   = [ 1 0 1 2   1 0 1 2   2 1 2 3   0 1 ];
states.multiplet = [ 1 1 1 1   2 2 2 2   3 3 3 3   4 4 ];
states.multipletNames = {'(s/d)1', 'p1', '(p/f)2', 'p0'};

if ~mix
    assign = { 0,1,1;   1,0,1;   1,1,1;   1,2,1; ...
               0,1,2;   1,0,2;   1,1,2;   1,2,2; ...
               0,2,1;   1,1,3;   1,2,3;   1,3,1; ...
               0,0,1;   1,1,6 };
    ref = [-0.048430, -0.009419];
else
    % Proposed by proposeMixAssignment.m at was almost right
    % Only wrong assigned: a(11)=jh2(6) -> jh2(5)
    % Energy order at the mix best fit ref (A=-0.0855, B=-0.0164), checked with
    % the channel weights. With mixing the quarkonium levels lie between
    % the hybrid ones, so some indices depend on A and B:
    %  - a(6)  = jh0(5): P1^{0+} (decoupled block). jh0(4) is the R0 3S
    %          quarkonium (94% Q) at this (A,B); it was jh0(4) at A=-0.0807
    %  - a(4)  = jh2(2) ((s/d)1 2-+, P1^{--}, mixing block) and
    %    a(8)  = jh2(3) (p1 2+-, decoupled block) are only 9 MeV apart and
    %          swap order around A=-0.09: the indices are the ones at the
    %          best fit (at A=-0.0807 they were jh2(3) and jh2(2))
    %  - a(12) = jh3(2) ((p/f)2 3+-, P1^{--}, mixing block) is 11 MeV from
    %          a decoupled level
    assign = { 0,1,6;   1,0,3;   1,1,3;   1,2,2; ...
               0,1,9;   1,0,5;   1,1,4;   1,2,3; ...
               0,2,6;   1,1,5;   1,2,5;   1,3,2; ...
               0,0,4;   1,1,9 };
    ref = [-0.085543, -0.016410];
end

if size(assign,1) ~= numel(t)
    error('hybridAssignments:size', 'assign has %d rows but there are %d data points.', ...
        size(assign,1), numel(t));
end
end
