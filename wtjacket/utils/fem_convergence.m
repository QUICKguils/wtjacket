function fem_convergence
% CONVERGENCE  Analyze the convergence of the FEM simulation.

% TODO:
% - make an opts to write in MAT file ?
% - analyse the NX modes, I don't know how to load them.
% - implement other ways to analyze the convergence, e.g. MAC w.r.t. the
%   NX modes ?

MLSol = gen_matlab_sol(1:5);
NXSol = load_nx_sol();

plot_freq_convergence(MLSol);
plot_freq_convergence(NXSol);

plot_mode_convergence(MLSol);
plot_mode_convergence(NXSol);
end

function MLSol = gen_matlab_sol(sdiv_set)
% GEN_MATLAB_SOL  Generate the solutions for a set of subdivisions.
%
% Argument:
%	sdiv_set (1xN int, default: 1:8) -- Set of desired number of subdivisions.

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));
% Results directory location.
res_file = (fullfile(file_dir, "../../res/"));

% Run the bare structure simulation once, if SOL data are not available.
if ~isfile(fullfile(res_file, "modeling_sol.mat"))
	modeling(1, '');
end
SOL = load(fullfile(res_file, "modeling_sol.mat"));

% Desired sdiv sequence.
if nargin == 0
	sdiv_set = 1:8;
end

% Gather the frequencies for the desired subdivisions.
freq_set = zeros(numel(sdiv_set), SOL.nbMode);
for i = 1:numel(sdiv_set)
	modeling(sdiv_set(i), 'w');
	SOL = load("res\modeling_sol.mat");
	freq_set(i, :) = SOL.freqs;
end

% Build return data structure.
MLSol.sdiv_set = sdiv_set;
MLSol.freq_set = freq_set;
MLSol.name     = 'Matlab';
end

function NXSol = load_nx_sol
% LOAD_NX_SOL  Load the solution found via the NX simulations.

% Sequence of number of subdivisions.
NXSol.sdiv_set = 1:8;

NXSol.freq_set = [ ...
	4.42401783E-01, 4.49817787E-01, 9.60608860E-01, 6.89354897E+00, 7.32969389E+00, 1.63825977E+01, 2.04311445E+01, 2.22578560E+01; ...
	4.42402461E-01, 4.49817194E-01, 9.60603890E-01, 6.87080562E+00, 7.30317814E+00, 1.63451574E+01, 1.97419262E+01, 2.14418085E+01; ...
	4.42402506E-01, 4.49817159E-01, 9.60603607E-01, 6.86900918E+00, 7.30109435E+00, 1.63422428E+01, 1.96804680E+01, 2.13663317E+01; ...
	4.42402515E-01, 4.49817152E-01, 9.60603550E-01, 6.86856939E+00, 7.30058720E+00, 1.63415557E+01, 1.96660325E+01, 2.13486387E+01; ...
	4.42402518E-01, 4.49817150E-01, 9.60603529E-01, 6.86839930E+00, 7.30039186E+00, 1.63412979E+01, 1.96606482E+01, 2.13420564E+01; ...
	4.42402520E-01, 4.49817149E-01, 9.60603520E-01, 6.86831582E+00, 7.30029625E+00, 1.63411740E+01, 1.96580727E+01, 2.13389142E+01; ...
	4.42402521E-01, 4.49817148E-01, 9.60603515E-01, 6.86826852E+00, 7.30024217E+00, 1.63411049E+01, 1.96566395E+01, 2.13371681E+01; ...
	4.42402521E-01, 4.49817148E-01, 9.60603513E-01, 6.86823904E+00, 7.30020852E+00, 1.63410622E+01, 1.96557576E+01, 2.13360949E+01];

NXSol.name = 'NX';
end

function plot_freq_convergence(SolSet)
% PLOT_FREQ_CONVERGENCE  Generate convergence graphs for the frequencies.
%
% Arguments:
%	SolSet (struct) -- Set of solution, with fields:
%	  sdiv_set (1xN int)    -- Set of number of subdivisions.
%	  freq_set (1xN double) -- Set of associated frequencies.
%	  name     (char)       -- Name of the solution.

% Build the solution name.
if SolSet.name ~= ""
	SolSet.name = " (" + SolSet.name + ")";
end

% Instantiate a figure object.
figure("WindowStyle", "docked");

% Plot the frequencies.
subplot(1, 2, 1);
plot(SolSet.sdiv_set, SolSet.freq_set);
title("Frequencies convergence" + SolSet.name);
xlabel("Number of sub-elements");
ylabel("Natural frequency (Hz)");
grid;

% Plot the residuals.
subplot(1, 2, 2);
% NOTE:
% this way of defining residuals is taken from the fluid mechanics
% course, section 3.3.
residuals = abs(diff(SolSet.freq_set)) ./ abs(SolSet.freq_set(2, :) - SolSet.freq_set(1, :));
semilogy(SolSet.sdiv_set(1:end-1), residuals);
title("Residuals convergence" + SolSet.name);
xlabel("Number of sub-elements");
ylabel("Residual");
grid;
end

function plot_mode_convergence(SolSet)
% TODO: implement this badboy
end
