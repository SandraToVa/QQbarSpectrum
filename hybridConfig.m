function cfg = hybridConfig(varargin)
% cfg = hybridConfig()                 returns the current settings
% cfg = hybridConfig('mix',true,...)   changes some settings (name-value pairs)
% cfg = hybridConfig('reset')          restores the defaults below
%
% Central place for the flags and parameters used by all the hybrid files
% GQQbarS0Jcal*.m and GQQbarS1Jcal*.m (their parameters2/3/4 read from here).
% The number of channels n is set automatically from the potential matrix,
% so switching mix on/off only requires changing it here (or calling
% hybridConfig('mix',true) at the top of a script).
%
% Settings persist until 'reset' or "clear functions"/"clear all".

persistent current
if isempty(current) || (nargin == 1 && strcmp(varargin{1}, 'reset'))
    current = defaults();
    if nargin == 1, cfg = current; return; end
end

for k = 1:2:nargin
    name = varargin{k};
    if ~isfield(current, name)
        error('hybridConfig:unknown', 'Unknown setting "%s".', name);
    end
    current.(name) = varargin{k+1};
end
cfg = current;
end


function d = defaults()
% Spectrum true value: the potentials give E, to obtain the spectrum we add
% 2mQ - Eval
d.Eg   = true;
d.Eval = 0.45;
% Mixing of hybrids with quarkonium
d.mix  = false;
% Hyperfine splitting (spin-1 hybrids)
d.hf     = true;
d.glamb1 = -0.1;   % glambda'   in GeV
d.glamb3 = 0.2;    % glambda''' in GeV
d.sig    = 0.21;   % string tension in GeV^2
d.pm     = -1;     % +- sign of Vsb
d.r0     = 3.96;   % GeV^-1
end
