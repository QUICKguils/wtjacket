function fem_convergence
% CONVERGENCE  Analyze the convergence of the FEM simulation.

% TODO:
% - analyse the NX modes, I don't know how to load them.
% - implement other ways to analyze the convergence, e.g. MAC w.r.t. the
%   NX modes ?

MLSol = gen_matlab_sol();
NXSol = load_nx_sol();

plot_freq_convergence(MLSol);
plot_freq_convergence(NXSol);

%plot_mass_convergence(MLSol);
%plot_mass_convergence(NXSol);

% plot_mode_convergence(MLSol);
% plot_mode_convergence(NXSol);
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
BS  = load(fullfile(res_file, "bare_struct.mat"));

% Desired sdiv sequence.
if nargin == 0
	sdiv_set = 1:8;
end

% Gather the frequencies, modes, masses for the desired subdivisions.
freq_set   = zeros(numel(sdiv_set), SOL.nbMode);
modes_set  = zeros(numel(sdiv_set), BS.nbDOF, SOL.nbMode);
masses_set = zeros(numel(sdiv_set), 2);

for i = 1:numel(sdiv_set)
	modeling(sdiv_set(i), '');
	SOL = load(fullfile(res_file, "modeling_sol.mat"));
	BS  = load(fullfile(res_file, "bare_struct.mat"));

	freq_set(i, :) = SOL.freqs;

	for d = 1:BS.nbDOF
		modes_set(i, d, :) = SOL.modes(d,:);
	end

	masses_set(i, 1) = SOL.mass_rbm;
	masses_set(i, 2) = SOL.mass_rbm/BS.mass;
end

% Build return data structure.
MLSol.sdiv_set   = sdiv_set;
MLSol.freq_set   = freq_set;
MLSol.modes_set  = modes_set;
MLSol.masses_set = masses_set;
MLSol.name = 'Matlab';
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

function plot_mass_convergence(SolSet)
% TODO mass convergence analysis

% Mass convergence
figure("WindowStyle", "docked");
subplot(1, 2, 1);
plot(sdiv_serie, SolSet.masses_set(:, 1));
title("RBM Total Mass convergence");
xlabel("Number of sub-elements");
ylabel("Total mass from RBM (kg)");
grid;

subplot(1, 2, 2);
plot(sdiv_serie, SolSet.masses_set(:, 2));
title("RBM Total Mass - Theoretical Mass Ratio convergence");
xlabel("Number of sub-elements");
ylabel("Mass ratio [-]");
grid;
end

function plot_mode_convergence(SolSet)
% Compute the relative differences between modes amplitudes
mode_diff = zeros(numel(sdiv_serie)-1, BS.nbDOF, SOL.nbMode);
norm_mode_diff = zeros(numel(sdiv_serie)-1, SOL.nbMode);

for i = 1:numel(sdiv_serie)-1
	mode_diff(i, :, :) = abs(modes_serie(i+1, :, :) -  modes_serie(i, :, :));

	for m = 1:SOL.nbMode
		norm_mode_diff(i,m) = norm(mode_diff(i, :, m));
	end
end

% MAC
% this section computes the MAC for each successive mode
MAC_arr = zeros(numel(sdiv_serie)-1, SOL.nbMode, SOL.nbMode);

for i = 1:numel(sdiv_serie)-1
	MAC_arr(i, :, :) = MAC(modes_serie(i, :, :), modes_serie(i+1, :, :), SOL);
end

% NOTE pretty certain MAC(mi, mi+1) should converge to eye(nbMode)
% if modes do converge ?

% Instantiate a figure object.
figure("WindowStyle", "docked");

% Plot of MAC convergence
for i = 1:numel(sdiv_serie)-1
	subplot(1, numel(sdiv_serie), i);
	image(reshape(MAC_arr(i,:,:), SOL.nbMode, SOL.nbMode) ,'CDataMapping','scaled')
	colorbar
end

% Mode convergence
subplot(1, 3, 3);
plot(sdiv_serie(2:end), norm_mode_diff);
title("Mode L2 norm convergence");
xlabel("Number of sub-elements");
ylabel("Relative difference mode norm");
grid;
end

function mac = MAC(m1, m2, SOL)
% MAC  Compute the Modal Assurance Criterion matrix.
%
% Arguments
%   m1   1 x nbMode x nbDOF
%   m2   1 x nbMode x nbDOF
%
% Returns
%   mac  nbMode x nbMode

% mac(m1, m2)[i][j] = ( (m1[i]' * m2[j]) / (|m1| * |m2|) )^2
mac = zeros(SOL.nbMode, SOL.nbMode);
for i = 1:SOL.nbMode
	for j = 1:SOL.nbMode

		norm_m1 = norm(m1(1, :, i));
		norm_m2 = norm(m2(1, :, j));

		mac(i,j) = (1 / (norm_m1 * norm_m2) * m1(1, :, i) * m2(1, :, j)').^2;
	end
end
end
