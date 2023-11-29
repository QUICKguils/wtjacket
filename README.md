# wtjacket

This repository holds all the code written for the project carried out as part
of the Theory of Vibration course (MECA0029-1), academic year 2023-2024.

## Basic usage

From the top level directory (where this README lies), just run the main project
file:
```matlab
run wtjacket\MECA0029_Group_3.m
```
The execution of this main function add all the project files to the matlab
path, for the current matlab session. Assuming that the top level directory path
is called `rootDir`, it is still possible to add the project files manually to
the matlab path:
```matlab
addpath(genpath(rootDir));
```

## Advanced usage

### Options of the main file

The default code execution parameters are stored in the
`wtjacket/util/load_defaults.m` file.
It is however possible to launch the main function overriding these default
values.
The docstring of the main function gives an exhaustive usage explanation.
In particular, it is possible to decide whether the relevant data that were
computed should be saved. If yes, a `res/` directory will be created at top
level, and will store these data in appropriate `.MAT` files.

For example:
```matlab
% Fetch the default code execution parameters.
RunArg = load_defaults();

% Override the defaults:
RunArg.sdiv   = 4;    % 4 subdivisions in the bare structure,
RunArg.method = 'n';  % Only use Newmark to compute the transient response.
RunARg.opts   = 's';  % Only save the computed data, don't plot them.

% Launch the main function.
MECA0029_Group_3(RunArg);
```

### Executing the main parts of the project independently

Apart from the main function `MECA0029_Group_3.m`, the `wtjacket/`
directory contains the functions that implement the three main parts of the
project, namely `modeling.m`, `transient.m` and `reduction.m`. They are the
"main subfunctions" of the project, and they can be run independently.

The usage of these main subfunctions is quite flexible. The user can refer to
the corresponding docstring for an exhaustive explanation.
For example:
```matlab
% Load the appropriated, previously generated data.
Stm        = load("res\statement.mat");
SdivStruct = load("res\subdivisedStructure.mat");
AlgSys     = load("res\algebraicSystem.mat");
FemSol     = load("res\femSolution.mat");

% Fetch and override the default code execution parameters:
RunArg = load_defaults();
RunArg.nMode      = 6;         % Number of computed first modes.
RunARg.opts       = 'p';       % Enable plots creation.
RunARg.nodeLabels = [14, 21];  % Compute solution for nodes 14 and 21.

% Run the code for the transient part, based on a modal superposition of the
% first six modes, for the nodes 14 and 21, and with plots drawing enabled.
transient(RunArg, Stm, SdivStruct, AlgSys, FemSol);
```

### Executing the analysis files

The functions that lie in the `wtjacket/analysis/` directory are also executable
as is. These functions implement a collection of utilities that allow to analyze
the results obtained from the main subfunctions.

For example, the convergence graphs of the natural frequencies that have been
computed can be obtained as follows:
```matlab
% Analyze the results of the modeling part, for a given range of subdivisions in
% the bare structure, and for the first eight computed modes.
modeling_analysis(1:8, 8);
```

## Project architecture

- `wtjacket/`:
  - `MECA_0029_Group_3.m` triggers all the code of the project.
  - `modeling.m` implements the first part of the project, namely the model of
    the wind turbine jacket, using 3D beam elements.
  - `transient.m` implements the second part of the project, namely the
    transient response due to a harmonic excitation.
  - `reduction.m` implements the third part of the project, namely the ...
  - `analysis/`: functions that analyze the results obtained from the main
    subfunctions.
  - `object/`: object-oriented representation of the vibration problem. In this
    directory are defined, among other things, classes that define nodes,
    elements, loads or concentrated masses.
  - `util/`: collection of utility functions that are used throughout the
    project.
- `res/`: contains the `.MAT` files that hold the relevant data that were
  computed by the main subfunctions. This directory is automatically created,
  if absent on the user's machine.
