function constants
% CONSTANTS  Constant quantities used thourghout the project.
%
% This function defines the constants used throughout the project, and
% write them in `constants.mat`.

%% Geometry

C.angle    = 3;                  % Angle of the main legs w.r.t. the Z-axis [°].
C.b_width  = 5;                  % Width at the base of the structure [m].
C.f_height = [0, 1, 9, 17, 25];  % Heights of the horizontal frames [m].
C.n_height = 80;                 % Height of the nacelle [m].


%% Save data into constants.mat

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Save data in constants.mat, which lies in the /res directory.
save(fullfile(file_dir, "../../res/constants.mat"), '-struct', "C");

end
