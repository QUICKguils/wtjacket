function MECA0029_Ernotte(varargin)
% MECA_0029_ERNOTTE  triggers all the code of the project.
%
% Arguments:
%	sdiv (int) -- Optional, default is 3.
%	  Number of subsivisions in the bare structure.
%	plt (char {'p', 'w'}) -- Optional, default is 'p'.
%	  Plotting options.
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.

%% Options setting

% Check number of inputs.
if numel(varargin) > 2
	error('At most 2 optional inputs are required');
end

% Set default value for optional inputs.
optargs = {3, 'p'};

% Overwrite default value of optional inputs.
optargs(1:numel(varargin)) = varargin;

% Place optional args in memorable variable names.
[sdiv, plt] = optargs{:};

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
bare_struct(plt);

%% Execute the code

% 1. Modeling of the structure.
modeling(sdiv, plt);

% 2. Transient response.
transient(plt);

% 3. Reduction methods.
reduction(plt);

end
