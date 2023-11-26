function modeling_analysis(varargin)
% MODELING_ANALYSIS  Analyze the results of the modeling part.
%
% Arguments:
%	sdivSet (1 x N int) -- Set of desired number of subdivisions (default: 1:8).
%	nMode   (int)       -- Number of first computed modes (default: 8).

% Set default value for optional inputs.
optargs = {1:8, 8};
% Overwrite default value of optional inputs.
optargs(1:numel(varargin)) = varargin;
% Place optional args in memorable variable names.
[sdivSet, nMode] = optargs{:};

MatlabSolution = compute_matlab_solution(sdivSet, nMode);
NxSolution = load_nx_solution(sdivSet, nMode);

plot_frequency_convergence(MatlabSolution);
plot_frequency_convergence(NxSolution);

% plot_mass_convergence(MLSol);
% plot_mass_convergence(NXSol);

% plot_mode_convergence(MLSol);
% plot_mode_convergence(NXSol);
end

function Solution = compute_matlab_solution(sdivSet, nMode)
% COMPUTE_MATLAB_SOLUTION  Compute the Matlab solutions.
%
% Arguments:
%	sdivSet (1 x N int) -- Set of desired number of subdivisions.
%	nMode   (int)       -- Number of first mode computed.
% Return:
%	Solution (struct) -- Set of solution, with fields:
%	  nMode   (int)          -- Number of first mode computed.
%	  sdivSet (1 x N int)    -- Set of number of subdivisions.
%	  freqSet (1 x N double) -- Set of associated frequencies.
%	  name    (char)         -- Name of the solution.

% Gather the frequencies, modes and masses for the desired subdivisions.
freqSet = zeros(numel(sdivSet), nMode);
% modeSet = zeros(numel(sdivSet), BS.nbDOF, nMode);
% massSet = zeros(numel(sdivSet), 2);

% Project statement data.
Stm = load_statement();

for i = 1:numel(sdivSet)
	[~, ~, ~, FemSol] = modeling(Stm, sdivSet(i), nMode, '');

	freqSet(i, :) = FemSol.frequencyHertz;

	% 	for d = 1:BS.nbDOF
	% 		modeSet(i, d, :) = SOL.modes(d,:);
	% 	end
	%
	% 	massSet(i, 1) = SOL.massFromRbm;
	% 	massSet(i, 2) = SOL.massFromRbm/BS.mass;
end

% Build return data structure.
Solution.sdivSet = sdivSet;
Solution.freqSet = freqSet;
% Solution.modeSet = modeSet;
% Solution.massSet = massSet;
Solution.nMode   = nMode;
Solution.name    = 'Matlab';
end

function Solution = load_nx_solution(sdivSet, nMode)
% LOAD_NX_SOLUTION  Load the solution found via the NX simulations.
%
% Arguments:
%	sdivSet (1xN int) -- Set of desired number of subdivisions.
%	nMode   (int)     -- Number of first mode computed.
% Return:
%	Solution (struct) -- Set of solution, with fields:
%	  nMode   (int)        -- Number of first mode computed.
%	  sdivSet (1xN int)    -- Set of number of subdivisions.
%	  freqSet (1xN double) -- Set of associated frequencies.
%	  name    (char)       -- Name of the solution.

% Computed on NX.
COMPUTED_FREQUENCIES = [ ...
	4.42401783E-01, 4.49817787E-01, 9.60608860E-01, 6.89354897E+00, 7.32969389E+00, 1.63825977E+01, 2.04311445E+01, 2.22578560E+01; ...
	4.42402461E-01, 4.49817194E-01, 9.60603890E-01, 6.87080562E+00, 7.30317814E+00, 1.63451574E+01, 1.97419262E+01, 2.14418085E+01; ...
	4.42402506E-01, 4.49817159E-01, 9.60603607E-01, 6.86900918E+00, 7.30109435E+00, 1.63422428E+01, 1.96804680E+01, 2.13663317E+01; ...
	4.42402515E-01, 4.49817152E-01, 9.60603550E-01, 6.86856939E+00, 7.30058720E+00, 1.63415557E+01, 1.96660325E+01, 2.13486387E+01; ...
	4.42402518E-01, 4.49817150E-01, 9.60603529E-01, 6.86839930E+00, 7.30039186E+00, 1.63412979E+01, 1.96606482E+01, 2.13420564E+01; ...
	4.42402520E-01, 4.49817149E-01, 9.60603520E-01, 6.86831582E+00, 7.30029625E+00, 1.63411740E+01, 1.96580727E+01, 2.13389142E+01; ...
	4.42402521E-01, 4.49817148E-01, 9.60603515E-01, 6.86826852E+00, 7.30024217E+00, 1.63411049E+01, 1.96566395E+01, 2.13371681E+01; ...
	4.42402521E-01, 4.49817148E-01, 9.60603513E-01, 6.86823904E+00, 7.30020852E+00, 1.63410622E+01, 1.96557576E+01, 2.13360949E+01];

% Extract desired frequencies.
Solution.freqSet = COMPUTED_FREQUENCIES(sdivSet, 1:nMode);

Solution.sdivSet = sdivSet;
Solution.nMode   = nMode;
Solution.name    = 'NX';
end

function plot_frequency_convergence(Solution)
% PLOT_FREQUENCY_CONVERGENCE  Generate convergence graphs for the frequencies.
%
% Arguments:
%	Solution (struct) -- Set of solution, with fields:
%	  nMode   (int)        -- Number of first mode computed.
%	  sdivSet (1xN int)    -- Set of number of subdivisions.
%	  freqSet (1xN double) -- Set of associated frequencies.
%	  name    (char)       -- Name of the solution.

% Build the solution name.
if Solution.name ~= ""
	solutionName = " (" + Solution.name + ")";
else
	solutionName = " (unknown solution) ";
end

% Instantiate a figure object.
figure("WindowStyle", "docked");

% Plot the frequencies.
subplot(1, 2, 1);
plot(Solution.sdivSet, Solution.freqSet);
title("Frequencies convergence" + Solution.name);
xlabel("Number of sub-elements");
ylabel("Natural frequency (Hz)");
grid;

% Plot the residuals.
subplot(1, 2, 2);
% NOTE:
% this way of defining residuals is taken
% from the fluid mechanics course, section 3.3.
residuals = abs(diff(Solution.freqSet)) ./ abs(Solution.freqSet(2, :) - Solution.freqSet(1, :));
semilogy(Solution.sdivSet(1:end-1), residuals);
title("Residuals convergence" + solutionName);
xlabel("Number of sub-elements");
ylabel("Residual");
grid;
end

function plot_mass_convergence(Solution)
% PLOT_MASS_CONVERGENCE  Plot the mass convergence.

% TODO mass convergence analysis

figure("WindowStyle", "docked");

subplot(1, 2, 1);
plot(Solution.sdivSet, Solution.massSet(:, 1));
title("RBM Total Mass convergence");
xlabel("Number of sub-elements");
ylabel("Total mass from RBM (kg)");
grid;

subplot(1, 2, 2);
plot(Solution.sdivSet, Solution.massSet(:, 2));
title("RBM Total Mass - Theoretical Mass Ratio convergence");
xlabel("Number of sub-elements");
ylabel("Mass ratio [-]");
grid;
end

function plot_mode_convergence()
% PLOT_MODE_CONVERGENCE  Compute the relative differences between modes amplitudes.

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
MAC_arr = zeros(numel(sdiv_serie)-1, SOL.nMode, SOL.nMode);

for i = 1:numel(sdiv_serie)-1
	MAC_arr(i, :, :) = compute_mac(modes_serie(i, :, :), modes_serie(i+1, :, :), SOL);
end

% NOTE pretty certain MAC(mi, mi+1) should converge to eye(nbMode)
% if modes do converge ?

% Instantiate a figure object.
figure("WindowStyle", "docked");

% Plot of MAC convergence
for i = 1:numel(sdiv_serie)-1
	subplot(1, numel(sdiv_serie), i);
	image(reshape(MAC_arr(i,:,:), SOL.nMode, SOL.nMode) ,'CDataMapping','scaled')
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

function mac = compute_mac(m1, m2, SOL)
% COMPUTE_MAC  Compute the Modal Assurance Criterion matrix.
%
% Arguments
%   m1   1 x nbMode x nbDOF
%   m2   1 x nbMode x nbDOF
%
% Returns
%   mac  nbMode x nbMode

% mac(m1, m2)[i][j] = ( (m1[i]' * m2[j]) / (|m1| * |m2|) )^2
mac = zeros(SOL.nMode, SOL.nMode);
for i = 1:SOL.nMode
	for j = 1:SOL.nMode

		norm_m1 = norm(m1(1, :, i));
		norm_m2 = norm(m2(1, :, j));

		mac(i,j) = (1 / (norm_m1 * norm_m2) * m1(1, :, i) * m2(1, :, j)').^2;
	end
end
end
