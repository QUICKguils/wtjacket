function MECA0029_Group_3(varargin)
% MECA0029_Group_3  Trigger all the code of the project.
%
% Arguments:
%	sdiv   (int)      -- Number of subsivisions in the bare structure (default: 3).
%	nMode  (int)      -- Number of computed first modes (default: 8).
%	method (1xN char) -- Methods used for the transient response (default: 'dan').
%	  'd' -> Mode [D]isplacement method.
%	  'a' -> Mode [A]cceleration method.
%	  'n' -> [N]ewmark (time integration).
%	opts   (1xN char) -- Options (default: 'ps').
%	  'p' -> Enable [P]lots creation.
%	  's' -> [S]ave generated data.

%% Options setting

% Set default value for optional inputs.
optargs = {3, 8, 'dan', 'ps'};

% Check number of inputs.
if nargin > numel(optargs)
	error('At most %u optional inputs are required', numel(optargs));
end

% Overwrite default value of optional inputs.
optargs(1:nargin) = varargin;

% Place optional args in memorable variable names.
[sdiv, nMode, method, opts] = optargs{:};

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
Cst = load_constants();

% 1. Modeling of the structure.
[BareStruct, SdivStruct, AlgSys, FemSol] = modeling(Cst, sdiv, nMode, opts);

% 2. Transient response.
[AlgSys, TransientSol] = transient(Cst, SdivStruct, AlgSys, FemSol, nMode, method, opts);

% 3. Reduction methods.
[frequencyHertz] = reduction(SdivStruct, AlgSys, nMode);

%% Save generated data

if contains(opts, 's')
	save(fullfile(resDirectory, "constants.mat"),           "-struct", "Cst");
	save(fullfile(resDirectory, "bareStructure.mat"),       "-struct", "BareStruct");
	save(fullfile(resDirectory, "subdivisedStructure.mat"), "-struct", "SdivStruct");
	save(fullfile(resDirectory, "algebraicSystem.mat"),     "-struct", "AlgSys");
	save(fullfile(resDirectory, "FemSolution.mat"),         "-struct", "FemSol");
	save(fullfile(resDirectory, "transientSolution.mat"),   "-struct", "TransientSol");
end

end
