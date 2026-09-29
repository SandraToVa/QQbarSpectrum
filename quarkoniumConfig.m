function cfg = quarkoniumConfig(varargin)
% cfg = quarkoniumConfig()                 returns the current settings
% cfg = quarkoniumConfig('mix',false,...)  changes some settings (name-value pairs)
% cfg = quarkoniumConfig('reset')          restores the defaults below
%
% Central place for the flags used by all the quarkonium files
% QQbarS0J*.m and QQbarS1J*.m (their parameters2/3 read from here).
% The number of channels n is set automatically from the potential matrix.
% Hybrid files use hybridConfig.m instead.
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
        error('quarkoniumConfig:unknown', 'Unknown setting "%s".', name);
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
% Mixing of quarkonium with hybrids
d.mix  = true;
end
