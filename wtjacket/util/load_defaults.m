function Default = load_defaults
% LOAD_DEFAULTS  Default execution parameters used throughout the source code.

% Number of subdivisions in the bare structure.
Default.sdiv = 3;

% Number of computed first modes.
Default.nMode = 8;

% Sample used for time evolutions, in seconds.
%
% NOTE:
% Mode superposition methods require a time step not larger than 0.002s,
% in order to obtain a coherent convergence. However, this default time
% sample gives satisfying results.
% See: analysis/transient_analysis.m
Default.tSet = 0:0.01:10;

% Label list of nodes to inspect.
Default.nodeLabels = [18, 22];

% Methods used to compute the transient response.
%   'd' -> Mode [D]isplacement method.
%   'a' -> Mode [A]cceleration method.
%   'n' -> [N]ewmark (time integration).
Default.method = 'dan';

% Number of first modes used in reductions.
Default.m = 3;

% Output options.
%   'p' -> Enable [P]lots creation.
%   's' -> [S]ave generated data.
Default.opts = 'ps';
end
