function MECA0029_Group_3(varargin)
% MECA0029_Group_3  triggers all the code of the project.
%
% Arguments:
%	sdiv  (int)      -- Number of subsivisions in the bare structure (default: 3).
%	nMode (int)      -- Number of computed first modes (default: 8).
%	opts  (1xN char) -- Options (default: 'ps').
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
%	  's' -> Save generated data.

%% Options setting

% Check number of inputs.
if numel(varargin) > 3
	error('At most 3 optional inputs are required');
end

% Set default value for optional inputs.
optargs = {3, 8, 'ps'};

% Overwrite default value of optional inputs.
optargs(1:numel(varargin)) = varargin;

% Place optional args in memorable variable names.
[sdiv, nMode, opts] = optargs{:};

%% Set program initial state

% Close previous plots.
close all;

% Find the root directory of the project.
rootDirectory = fullfile(fileparts(mfilename('fullpath')), "..");

% Create the untracked results directory, if absent.
resDirectory = fullfile(rootDirectory, "/res");
if ~isfolder(resDirectory)
	mkdir(resDirectory);
end

% Add resursively sub-directories in the Matlab path.
addpath(genpath(fullfile(rootDirectory, "wtjacket")));

%% Execute the code

% 0. Load the constants.
C = load_constants();

% 1. Modeling of the structure.
[BS, SS, KM, SOL] = modeling(C, sdiv, nMode, opts);

% 2. Transient response.
% transient(C, KM, SOL, opts);

% 3. Reduction methods.
reduction(opts);

%% Save generated data

if contains(opts, 's')
	save(fullfile(resDirectory, "constants.mat"),           "-struct", "C");
	save(fullfile(resDirectory, "bareStructure.mat"),       "-struct", "BS");
	save(fullfile(resDirectory, "subdivisedStructure.mat"), "-struct", "SS");
	save(fullfile(resDirectory, "globalMatrices.mat"),      "-struct", "KM");
	save(fullfile(resDirectory, "modelingSolution.mat"),    "-struct", "SOL");
end

end
