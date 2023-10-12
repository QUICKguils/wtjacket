function constants
% CONSTANTS  Constant quantities used thourghout the project.
%
% This function defines the constants used throughout the project, and
% write them in `constants.mat`.

%% Geometry

C.leg_angle       = 3;                  % Angle of the main legs w.r.t. the Z-axis [°].
C.base_width      = 5;                  % Width at the base of the structure [m].
C.frame_height    = [0, 1, 9, 17, 25];  % Heights of the horizontal frames [m].
C.nacelle_height  = 80;                 % Height of the nacelle [m].

%% Transient response

C.impact_freq       = 1;         % Frequency of the resulting force on the jacket [Hz].
C.tail_weight       = 1e3;       % Weight of one whale tail [kg].
C.impact_speed      = 25 / 3.6;  % Velocity of the tail during the impact [m/s].
C.impact_duration   = 0.05;      % Duration of the impact [s].
C.momentum_transfer = 0.85;      % Proportion of the momentum transferred to the jacket [-].
C.force_direction   = 45;        % Direction of the impact force, w.r.t. the X-axis.
C.damping_ratio     = 5e-3;      % Damping ratio of the first two modes [-].

%% Save data into constants.mat

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Save data in constants.mat, which lies in the /res directory.
save(fullfile(file_dir, "../../res/constants.mat"), '-struct', "C");

end
