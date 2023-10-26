function Cst = load_constants
% LOAD_CONSTANTS  Constants quantities used thourghout the project.
%
% This function returns a structure that contains the constant
% quantities used throughout the project.

%% Geometry

Cst.LEG_ANGLE      = 3;                  % Angle of the main legs w.r.t. the Z-axis [°].
Cst.BASE_WIDTH     = 5;                  % Width at the base of the structure [m].
Cst.FRAME_HEIGHT   = [0, 1, 9, 17, 25];  % Heights of the horizontal frames [m].
Cst.NACELLE_HEIGHT = 80;                 % Height of the nacelle [m].

%% Transient response

Cst.LOAD_FREQUENCY_HERTZ = 1;         % Frequency of the resulting force on the jacket [Hz].
Cst.TAIL_MASS            = 1e3;       % Mass of one whale tail [kg].
Cst.TAIL_SPEED           = 25 / 3.6;  % Velocity of the tail during the impact [m/s].
Cst.IMPACT_DURATION      = 0.05;      % Duration of the impact [s].
Cst.MOMENTUM_TRANSFER    = 0.85;      % Proportion of the momentum transferred to the jacket [-].
Cst.FORCE_DIRECTION      = 45;        % Direction of the impact force, w.r.t. the X-axis [°].
Cst.INITIAL_CONDITIONS   = [0, 0];    % Initial conditions of the transient response [m, m/s].
Cst.DAMPING_RATIO        = 5e-3;      % Damping ratio of the first two modes [-].

end
