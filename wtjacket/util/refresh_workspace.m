% REFRESH_WORKSPACE  Bring computed data in global workspace.

% Find the results directory
resDirectory = fullfile(fileparts(mfilename('fullpath')), "../../res/");

RunArg       = load(fullfile(resDirectory, "runArguments.mat"));
Stm          = load(fullfile(resDirectory, "statement.mat"));
BareStruct   = load(fullfile(resDirectory, "bareStructure.mat"));
SdivStruct   = load(fullfile(resDirectory, "subdivisedStructure.mat"));
AlgSys       = load(fullfile(resDirectory, "algebraicSystem.mat"));
FemSol       = load(fullfile(resDirectory, "femSolution.mat"));
TransientSol = load(fullfile(resDirectory, "transientSolution.mat"));
ReductionSol = load(fullfile(resDirectory, "reductionSolution.mat"));

clear resDirectory;
