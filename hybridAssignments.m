function [assign, t, e] = hybridAssignments(mix)
% [assign, t, e] = hybridAssignments(mix)
% Lattice data used in QQbarSpectrum_Optimize and the computed state that
% each data point is compared with. Row k of assign corresponds to t(k), e(k).
%
% assign{k,:} = {spin, Jcal, index}
%   spin  : 0 -> GQQbarS0Jcal<Jcal>.m,  1 -> GQQbarS1Jcal<Jcal>.m
%   index : position of the state in energy order (1 = lowest) within that
%           file, at the best-fit (A,B). The optimizer checks this
%           self-consistently: if the order at the minimum differs from the
%           one used during the fit, it refits with the new order.
%
% The indices change with mix (quarkonium states appear in the spectrum),
% so there is one table per case. proposeMixAssignment.m suggests the
% mix=true table from the mix=false one (overlap of wave functions).

% Lattice data (order: {4 (s/d)1 states}, {4 p_1 states}, {4 (p/f)2 states},
% {2 p_0 states})
t = [4.0296 3.8976 3.9286 4.0746 4.1436 4.1106 4.1116 4.1756 ...
     4.2306 4.1786 4.2386 4.2516 4.4396 4.5136];
e = [0.0176 0.0186 0.0236 0.0216 0.0256 0.0276 0.0236 0.0186 ...
     0.0326 0.0276 0.0266 0.0346 0.0466 0.0536];

if ~mix
    assign = { 0,1,1;   1,0,1;   1,1,1;   1,2,1; ...
               0,1,2;   1,0,2;   1,1,2;   1,2,2; ...
               0,2,1;   1,1,3;   1,2,4;   1,3,1; ...
               0,0,1;   1,1,6 };
else
    % Proposed by proposeMixAssignment.m at A=-0.0548, B=0.0038 (overlaps
    % 0.95-1.00 with the mix=false states). Row 13: the mix=false channel
    % of GQQbarS0Jcal0 is not present with mixing; state 4 is the one that
    % is 97% Sigma_u hybrid. Check it if the mix=true results look odd.
    assign = { 0,1,6;   1,0,3;   1,1,3;   1,2,2; ...
               0,1,9;   1,0,4;   1,1,4;   1,2,3; ...
               0,2,6;   1,1,5;   1,2,6;   1,3,2; ...
               0,0,4;   1,1,9 };
end

if size(assign,1) ~= numel(t)
    error('hybridAssignments:size', 'assign has %d rows but there are %d data points.', ...
        size(assign,1), numel(t));
end
end
