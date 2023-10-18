function C = load_constants
% LOAD_CONSTANTS  Constants quantities used thourghout the project.
%
% This function returns a structure that contains the constant
% quantities used throughout the project.

%% Geometry

C.LEG_ANGLE      = 3;                  % Angle of the main legs w.r.t. the Z-axis [°].
C.BASE_WIDTH     = 5;                  % Width at the base of the structure [m].
C.FRAME_HEIGHT   = [0, 1, 9, 17, 25];  % Heights of the horizontal frames [m].
C.NACELLE_HEIGHT = 80;                 % Height of the nacelle [m].

%% Transient response

C.LOAD_FREQUENCY    = 1;         % Frequency of the resulting force on the jacket [Hz].
C.TAIL_MASS         = 1e3;       % Mass of one whale tail [kg].
C.TAIL_SPEED        = 25 / 3.6;  % Velocity of the tail during the impact [m/s].
C.IMPACT_DURATION   = 0.05;      % Duration of the impact [s].
C.MOMENTUM_TRANSFER = 0.85;      % Proportion of the momentum transferred to the jacket [-].
C.FORCE_DIRECTION   = 45;        % Direction of the impact force, w.r.t. the X-axis [°].
C.DAMPING_RATIO     = 5e-3;      % Damping ratio of the first two modes [-].

end
