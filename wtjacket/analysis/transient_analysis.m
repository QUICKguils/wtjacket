function transient_analysis(varargin)
	% TRANSIENT_ANALYSIS  Analyze the results of the transient part.
	%
	% Arguments:
	%	focus (1xN char) -- Analysis to focus on (default: 'lsn').
	%	  'l' -> [L]oad spatial distribution.
	%	  's' -> [S]teady-state response.
	%	  'c' -> [C]onvergence of the mode superposition methods.
	%	  'n' -> [N]ewmark solution.

	% TODO:
	% - better management of sdiv and nMode across functions.

	% Set default value for optional inputs.
	optargs = {'lcsn', 'da'};
	% Overwrite default value of optional inputs.
	optargs(1:numel(varargin)) = varargin;
	% Place optional args in memorable variable names.
	[focus, cMethod] = optargs{:};

	% TODO: remove that when debugged.
	close all;

	% Compute the FEM solution.
	Cst = load_constants();
	[~, SdivStruct, AlgSys, FemSol] = modeling(Cst, 3, 8, '');

	% Wrapper for the transient response computations.
	this_transient = @(nMode, method, opts) transient(Cst, SdivStruct, AlgSys, FemSol, nMode, method, opts);

	% Analyze the transient code.
	if contains(focus, 'l'); analyze_load_distribution(FemSol, this_transient); end
	if contains(focus, 's'); analyze_steady_state(SdivStruct.nodeList, FemSol, this_transient); end
	if contains(focus, 'c'); analyze_convergence(); end
	if contains(focus, 'n'); analyze_NewmarkSol(SdivStruct.nodeList, this_transient); end
end

%% 1. Analyze the load spatial distribution

function analyze_load_distribution(FemSol, this_transient)
	% ANALYZE_LOAD_DISTRIBUTION  Spatial factor of the M-norm of the displacements.
	%
	% Arguments:
	%	FemSol         (struct) -- Solution of the FEM simulation.
	%	this_transient (handle) -- Compute the desired transient solution.

	% Get the modal superposition parameters.
	[~, TransientSol] = this_transient(8, '', '');
	% Local aliases, for readability.
	g  = TransientSol.DiscreteLoad.spatial;
	mu = TransientSol.ModalSup.mu;
	x  = FemSol.mode;

	% Spatial factor: related to max. of modal participation factors phi.
	qSpatialNorm = abs(x'*g) ./ sqrt(mu);

	figure("WindowStyle", "docked");
	stem(qSpatialNorm);
	set(gca,'yscal','log');
	grid;
	title([ ...
	"Distribution of the ||q||M spatial factors", ...
	"across the first " + num2str(TransientSol.ModalSup.nMode) + " modes"]);
	xlabel("Associated mode number");
	ylabel("Displacement spatial M-norms");
end

%% 2. Analyze the steady-state displacement response

function analyze_steady_state(nodeList, FemSol, this_transient)
	% Plot the steady-state displacement amplitude.
	%
	% This function uses the spectral expansion of the admittance matrix to compute
	% the amplitude of the displacement.
	% It thus assumes a harmonic excitation of frequency `w`.
	%
	% Arguments:
	%	nodeList       {1xN Node} -- Cell list of nodes.
	%	FemSol         (struct)   -- Solution of the FEM simulation.
	%	this_transient (handle)   -- Compute the desired transient solution.

	% TODO: make this function more flexible, not only (18, 18)

	% Get the modal superposition parameters.
	[AlgSys, TransientSol] = this_transient(8, '', '');
	% Local aliases, for readability.
	g  = TransientSol.DiscreteLoad.spatial;
	mu = TransientSol.ModalSup.mu;

	function H = build_H(AlgSys, FemSol, mu, wSet)
		% Build the admittance matrix H, for the given frequencies.
		%
		% Argument:
		%	wSet (1xN double) -- Harmonic excitation frequencies [rad/s].
		% Return:
		%	H (nDofFreexnDofFree complex) -- Admittance matrix.

		% Preallocate admittance matrix.
		H = zeros(AlgSys.nDofFree, AlgSys.nDofFree, numel(wSet));
		% Local aliases, for readability.
		x   = FemSol.mode;
		eps = AlgSys.eps;
		w0  = FemSol.frequencyRad;
		% Precompute the mode dyadic products.
		xTx = zeros(AlgSys.nDofFree, AlgSys.nDofFree, FemSol.nMode);
		for s = 1:FemSol.nMode
			xTx(:, :, s) = x(:, s) * x(:, s)';
		end
		% Build H through spectral expansion.
		for w = 1:numel(wSet)
			for s = 1:FemSol.nMode
				H(:, :, w) = H(:, :, w) + xTx(:, :, s) ./ (mu(s) * (w0(s)^2 - wSet(w)^2 + 2i*eps(s)*w0(s)*wSet(w)));
			end
		end
	end

	function qSteadyAmpl = compute_ss_amplitude(AlgSys, FemSol, mu, wSet)
		% COMPUTE_SS_AMPLITUDE  Compute the steady state response amplitude.
		%
		% This function computes the amplitude of the steady-state displacement,
		% for the given harmonic load frequency.

		HSet = build_H(AlgSys, FemSol, mu, wSet);
		qSteadyAmpl = zeros(AlgSys.nDofFree, numel(wSet));
		for w = 1:numel(wSet)
			qSteadyAmpl(:, w) = abs(HSet(:, :, w) * g);
		end
	end

	freqHz = 0:0.05:20;
	freqRad = freqHz * 2*pi;

	qSteadyAmpl = compute_steady_state_amplitude(AlgSys, FemSol, mu, freqRad);

	x_proj = qSteadyAmpl(nodeList{18}.dof(1), :) * cosd(45) + qSteadyAmpl(nodeList{18}.dof(2), :) * sind(45);

	figure("WindowStyle", "docked");
	semilogy(freqHz, x_proj);
	grid;
	xlabel("Frequency (Hz)");
	ylabel("Displacement (node 18, dir of load)");
	title("Amplitude of the steady harmonic load");
	% -> pic vers 4e et 5e comme attendu. Mais en fait plus gros pic encore vers
	% 0.44Hz, càd vers la première fréq. Si on essaye de changer la fréquence
	% d'excitation de 1hz à celle-ci, on voit qu'en effet la réponse transitoire
	% s'amplifie fort au cours du temps, jusqu'à atteindre une valeur fort
	% grande.

	xAt1 = compute_steady_state_amplitude(AlgSys, FemSol, mu, 2*pi);
	disp(xAt1(nodeList{18}.dof(1), :)*cosd(45) + xAt1(nodeList{18}.dof(2), :)*sind(45));
end

%% 3. Analyze the convergence of the mode superposition methods

function analyze_convergence(nModeSet, method)
	% ANALYZE_CONVERGENCE  Convergence of the mode superposition methods.
	%
	% Arguments:

	if contains(method, 'd'); dConv = this_method_convergence('d', nModeSet); end
	if contains(method, 'a'); aConv = this_method_convergence('a', nModeSet); end

	function this_method_convergence(method, nModeSet)
		% THIS_METHOD_CONVERVENCE  Convergence for the given method.
		nnMode = numel(nModeSet);
		for inMode = 1:nnMode
			compute_norm(method, nModeSet(inMode));
		end
	end

	function compute_norm(method, nMode)
		% COMPUTE_NORM  Compute the M-norm of the displacements.

		% Get the solution of the transient problem,
		% for the mode displacement and the mode acceleration methods.
		[AlgSys, TransientSol] = this_transient(nMode, method, '');
		% Local aliases, for readability.
		M = AlgSys.M_free;
		q = TransientSol.ModeDisplacement.q;

		% Compute the M-norm.
		qNorm = sqrt(q' * M * q);
		disp(qNorm);
	end

	function plot_convergence()

	end
end

%% 4. Analyze the Newmark solution

function analyze_NewmarkSol(nodeList, this_transient)
	% ANALYZE_NEWMARK  Analyze the Newmark solution.
%
% Arguments:
%	nodeList       {1xN Node} -- Cell list of nodes.
%	this_transient (handle)   -- Compute the desired transient solution.

	% WARN: the fft strongly depend on the chosen time interval.
% make sure to integrate enough 1st mode periods.

	% Get the solution of the transient problem,
% for the Newmark method.
	[~, TransientSol] = this_transient(8, 'n', '');
	% Local aliases, for readability.
	timeStep      = TransientSol.TimeParams.steps(1);  % WARN: not robust
	loadDirection = TransientSol.DiscreteLoad.direction;
	nMode         = TransientSol.ModalSup.nMode;
	q             = TransientSol.Newmark.q;

	inspectNodeLabels = [18, 22];

	figure("WindowStyle", "docked");

	nNode = numel(inspectNodeLabels);
	for iNode = 1:nNode
		% Displacement along the load direction.
		qX = q(nodeList{inspectNodeLabels(iNode)}.dof(1), :) * loadDirection(1);
		qY = q(nodeList{inspectNodeLabels(iNode)}.dof(2), :) * loadDirection(2);
		qZ = q(nodeList{inspectNodeLabels(iNode)}.dof(3), :) * loadDirection(3);
		qDir = qX + qY + qZ;

		% FFT of the displacement.
		y = fft(qDir);
		% Correct the scaling and shifting.
		fs = 1/timeStep;
		n = length(qDir);
		fshift = (-n/2:n/2-1) * (fs/n);
		yshift = fftshift(y);

		% Plot the FFT.
		subplot(nNode, 1, iNode);
		semilogy(fshift(ceil(end/2:end)), abs(yshift(ceil(end/2:end))));
		title('FFT of the Newmark''s transient response', ...
		['(node: ', num2str(inspectNodeLabels(iNode)), ...
		', order: ', num2str(nMode), ')']);
		xlabel("Frequency (Hz)");
		ylabel("FFT");
		xlim([0, 30]);
		grid;
	end
end
