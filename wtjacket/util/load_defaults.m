function Def = load_defaults
% LOAD_DEFAULTS  Default values used throughout the source code.

% Number of subdivisions in the bare structure.
Def.SDIV = 3;

% Number of computed first modes.
Def.N_MODE = 8;

% Sample used for time evolutions, in seconds.
Def.T_SET = 0:0.01:10;

% Methods used to compute the transient response.
%   'd' -> Mode [D]isplacement method.
%   'a' -> Mode [A]cceleration method.
%   'n' -> [N]ewmark (time integration).
Def.METHOD = 'dan';

% Output options.
%   'p' -> Enable [P]lots creation.
%   's' -> [S]ave generated data.
Def.OPTS   = 'ps';
end
