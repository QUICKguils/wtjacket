function Stm = load_statement
% LOAD_STATEMENT  Information provided by the project statement.
%
% This function returns a structure that contains the data provided by
% the statement, which are used throughout the project.

%% Geometry

Stm.LEG_ANGLE      = 3;                  % Angle of the main legs w.r.t. the Z-axis [°].
Stm.BASE_WIDTH     = 5;                  % Width at the base of the structure [m].
Stm.FRAME_HEIGHT   = [0, 1, 9, 17, 25];  % Heights of the horizontal frames [m].
Stm.NACELLE_HEIGHT = 80;                 % Height of the nacelle [m].

%% Transient response

Stm.LOAD_FREQUENCY_HERTZ = 1;         % Frequency of the resulting force on the jacket [Hz].
Stm.TAIL_MASS            = 1e3;       % Mass of one whale tail [kg].
Stm.TAIL_SPEED           = 25 / 3.6;  % Velocity of the tail during the impact [m/s].
Stm.IMPACT_DURATION      = 0.05;      % Duration of the impact [s].
Stm.MOMENTUM_TRANSFER    = 0.85;      % Proportion of the momentum transferred to the jacket [-].
Stm.FORCE_DIRECTION      = 45;        % Direction of the impact force, w.r.t. the X-axis [°].
Stm.INITIAL_CONDITIONS   = [0, 0];    % Initial conditions of the transient response [m, m/s].
Stm.DAMPING_RATIO        = 5e-3;      % Damping ratio of the first two modes [-].

end
