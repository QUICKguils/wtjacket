function constants
% CONSTANTS  Constant quantities used thourghout the project.
%
% This function defines the constants used throughout the project, and
% write them in `constants.mat`.

%% Geometry

C.leg_angle       = 3;                  % Angle of the main legs w.r.t. the Z-axis [°].
C.base_width      = 5;                  % Width at the base of the structure [m].
C.frame_height    = [0, 1, 9, 17, 25];  % Heights of the horizontal frames [m].
C.nacelle_height  = 80;                 % Height   of the nacelle [m].
C.nacelle_mass    = 200e3;              % Mass     of the nacelle [kg].
C.nacelle_inertia = 24e6;               % Inertias of the nacelle [kg*m²].

%% Transient response

C.ship_mass       = 50e6;  % Mass of the ship [kg].
C.ship_speed      = 0.25;  % Speed of the ship [m/s].
C.impact_duration = 1;     % Duration of the impact [s].
C.rem_momentum    = 0.8;   % Fraction of the remaining vs. initial momentum [-].
C.force_angle     = 45;    % Angle of the force w.r.t. the X-axis [°].
C.eps_s           = 0.03;  % Modal damping coefficient, for the eight first modes [-].

%% Save data into constants.mat

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Save data in constants.mat, which lies in the /res directory.
save(fullfile(file_dir, "../../res/constants.mat"), '-struct', "C");

end
