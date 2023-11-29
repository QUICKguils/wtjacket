function MECA0029_Group_3(RunArg)
% MECA0029_Group_3  Trigger all the code of the project.
%
% Argument:
%	RunArg (struct) -- Optional code execution parameters, with fields:
%	  sdiv   (int)        -- Number of subsivisions in the bare structure.
%	  nMode  (int)        -- Number of computed first modes.
%	  tSet   (1xN double) -- Time sample used for time evolutions.
%	  method (1xN char)   -- Methods used for the transient response.
%	    'd' -> Mode [D]isplacement method.
%	    'a' -> Mode [A]cceleration method.
%	    'n' -> [N]ewmark (time integration).
%	  opts   (1xN char) -- Output options.
%	    'p' -> Enable [P]lots creation.
%	    's' -> [S]ave generated data.
%
% The default values unsed to run this function
% are stored in util/load_defaults.m

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

%% Options setting

% Fetch the defaults execution parameters.
Default = load_defaults();

% Overwrite these defaults with user input.
switch nargin
	case 0
		RunArg = Default;
	case 1
		for fn = fieldnames(Default)'
			if ~isfield(RunArg, fn)
				RunArg.(fn{:}) = Default.(fn{:});
			end
		end
	otherwise
		error("Wrong number of input parameters.");
end

%% Execute the code

% 0. Load the project statement data.
Stm = load_statement();

% 1. Modeling of the structure.
[BareStruct, SdivStruct, AlgSys, FemSol] = modeling(RunArg, Stm);

% 2. Transient response.
[AlgSys, TransientSol] = transient(RunArg, Stm, SdivStruct, AlgSys, FemSol);

% 3. Reduction methods.
[GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol, CBReducedAlgSys, CBReducedFemSol, ReducedNewmarkSol] = reduction(Cst, SdivStruct, AlgSys, nMode, 3, opts);

%% Save generated data

if contains(RunArg.opts, 's')
	save(fullfile(resDirectory, "runArguments.mat"),        "-struct", "RunArg");
	save(fullfile(resDirectory, "statement.mat"),           "-struct", "Stm");
	save(fullfile(resDirectory, "bareStructure.mat"),       "-struct", "BareStruct");
	save(fullfile(resDirectory, "subdivisedStructure.mat"), "-struct", "SdivStruct");
	save(fullfile(resDirectory, "algebraicSystem.mat"),     "-struct", "AlgSys");
	save(fullfile(resDirectory, "femSolution.mat"),         "-struct", "FemSol");
	save(fullfile(resDirectory, "transientSolution.mat"),   "-struct", "TransientSol");
	save(fullfile(resDirectory, "GIReducedFemSol.mat"),		"-struct", "GIReducedFemSol");
	save(fullfile(resDirectory, "CBReducedFemSol.mat"),     "-struct", "CBReducedFemSol");
	save(fullfile(resDirectory, "ReducedNewmarkSol.mat"),   "-struct", "ReducedNewmarkSol");
end

end
