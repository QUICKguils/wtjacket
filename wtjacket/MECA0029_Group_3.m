function MECA0029_Group_3(varargin)
% MECA0029_Group_3  triggers all the code of the project.
%
% Arguments:
%	sdiv (int) -- Optional, default is 3.
%	  Number of subsivisions in the bare structure.
%	plt (char {'p', 'w'}) -- Optional, default is 'p'.
%	  Plotting options.
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.

% TODO:
% - will the write option eventually be written ?

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

% Create the /res untracked file, if absent.
res_dir = fullfile(root_dir, "/res");
if ~isfolder(res_dir)
	mkdir(res_dir);
end

% Add resursively sub-directories in the Matlab path.
addpath(genpath(fullfile(root_dir, "wtjacket")));

%% Execute the code

% 1. Modeling of the structure.
modeling(sdiv, plt);

% 2. Transient response.
transient(plt);

% 3. Reduction methods.
reduction(plt);

end
