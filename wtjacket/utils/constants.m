function constants
% CONSTANTS  Constant quantities used thourghout the project.
%
% This function defines the constants used throughout the project, and
% write them in `constants.mat`.

%%

C.rho = 7800;


%% Save data into constants.mat

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Save data in constants.mat, which lies in the root directory.
save(fullfile(file_dir, "../../res/constants.mat"), '-struct', "C");

end