function MECA0029_Group_3(varargin)
% MECA0029_Group_3  triggers all the code of the project.
%
% Arguments:
%	sdiv (int, default: 3)        -- Number of subsivisions in the bare structure.
%	opts (1xN char, default: 'p') -- Options.
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
%	  'w' -> Enable warnings.

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
[sdiv, opts] = optargs{:};

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

% Reset warnings state.
warning('on');

%% Execute the code

% 1. Modeling of the structure.
modeling(sdiv, opts);

% 2. Transient response.
transient(opts);

% 3. Reduction methods.
reduction(opts);

end