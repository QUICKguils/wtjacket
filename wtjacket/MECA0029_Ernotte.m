function MECA0029_Ernotte(opts)
% MECA_0029_ERNOTTE  triggers all the code of the project.
%
% Argument:
%	opts: char {'p', 'w'}, optional. Default is 'p'.
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.

close all;

%% Options setting

% Option defaults: generate the plots.
if ~nargin
	opts = 'p';
end

%% Set path and global MAT files

% Find the root directory of the project.
root_dir = fullfile(fileparts(mfilename('fullpath')), "..");

% % Add resursively sub-directories in the Matlab path.
addpath(genpath(fullfile(root_dir, "wtjacket")));

% Initialize MAT file.
constants();
bare_struct(opts);

%% Execute the code

% Load project constants.
C  = load(fullfile(root_dir, "res/constants.mat"));
BS = load(fullfile(root_dir, "res/bare_struct.mat"));

% 1. Modeling of the structure.
sdiv = 3;
modeling(sdiv, opts);

% 2. Transient response.

% 3. Reduction methods.

end