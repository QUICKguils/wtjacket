function MECA0029_Ernotte(sdiv, opts)
% MECA_0029_ERNOTTE  triggers all the code of the project.
%
% Arguments:
%	sdiv: int, optional. Default is 3.
%	  Number of subsivisions in the bare structure.
%	opts: char {'p', 'w'}, optional. Default is 'p'.
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.

%% Options setting

if nargin == 0
	sdiv = 3;
	opts = 'p';
elseif nargin == 1
	opts = 'p';
end

%% Set program initial state

% Find the root directory of the project.
root_dir = fullfile(fileparts(mfilename('fullpath')), "..");

% % Add resursively sub-directories in the Matlab path.
addpath(genpath(fullfile(root_dir, "wtjacket")));

% Reset class internal states, close previous plots.
clear Node Elem
close all;

% Initialize MAT file.
constants();
bare_struct(opts);

%% Execute the code

% Load project constants.
C  = load(fullfile(root_dir, "res/constants.mat"));
BS = load(fullfile(root_dir, "res/bare_struct.mat"));

% 1. Modeling of the structure.
modeling(sdiv, opts);

% 2. Transient response.

% 3. Reduction methods.

end