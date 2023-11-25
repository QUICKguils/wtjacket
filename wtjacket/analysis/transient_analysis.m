function transient_analysis(varargin)
	% TRANSIENT_ANALYSIS  Analyze the results of the transient part.
	%
	% Arguments:
	%	focus (1xN char) -- Analysis to focus on (default: 'lscn').
	%	  'l' -> [L]oad spatial distribution.
	%	  's' -> [S]teady-state response.
	%	  'c' -> [C]onvergence of the mode superposition methods.
	%	  'n' -> [N]ewmark solution.

	% TODO:
	% - better management of sdiv and nMode across functions.
	% - Maybe a util function for getting the projection of q for a given
	%   direction and from a given node.

	% Set default value for optional inputs.
	optargs = {'lscn'};
	% Overwrite default value of optional inputs.
	optargs(1:numel(varargin)) = varargin;
	% Place optional args in memorable variable names.
	focus = optargs{:};

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
	if contains(focus, 'c'); analyze_convergence(this_transient, SdivStruct.nodeList, 18, 1:8); end
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

function analyze_steady_state(nodeList, FemSol, this_transient, inspectNodeLabels)
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
	loadDirection = TransientSol.DiscreteLoad.direction;
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

	function qSteadyAmpl = compute_ss_amplitude(AlgSys, FemSol, mu, g, freqRadSet)
		% COMPUTE_SS_AMPLITUDE  Compute the steady state response amplitude.
		%
		% This function computes the amplitude of the steady-state displacement,
		% for the given harmonic load frequency.

		HSet = build_H(AlgSys, FemSol, mu, freqRadSet);
		qSteadyAmpl = zeros(AlgSys.nDofFree, numel(freqRadSet));
		for w = 1:numel(freqRadSet)
			qSteadyAmpl(:, w) = abs(HSet(:, :, w) * g);
		end
	end

	freqHzSet = 0:0.05:20;
	freqRadSet = freqHzSet * 2*pi;

	qSteadyAmpl = compute_ss_amplitude(H, g, freqRadSet);

	qAmplX = Method.q(nodeList{inspectNodeLabels(iNode)}.dof(1), :) * loadDirection(1);
	qAmplY = Method.q(nodeList{inspectNodeLabels(iNode)}.dof(2), :) * loadDirection(2);
	qAmplZ = Method.q(nodeList{inspectNodeLabels(iNode)}.dof(3), :) * loadDirection(3);

	qAmplDir = qAmplX + qAmplY + qAmplZ;


	x_proj = qSteadyAmpl(nodeList{18}.dof(1), :) * cosd(45) + qSteadyAmpl(nodeList{18}.dof(2), :) * sind(45);

	figure("WindowStyle", "docked");

	semilogy(freqHzSet, x_proj);
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

function analyze_convergence(this_transient, nodeList, nodeLabel, nModeSet)
	% ANALYZE_CONVERGENCE  Convergence of the mode superposition methods.
	%
	% Arguments:

	% TODO:
	% - maybe diffs the solution to a constant ref. For example the newmark sol.
	% - understant what the hell is going on with these convergences.

	nnMode = numel(nModeSet);

	% Get solution parameters at disposal.
	[~, TransientSol] = this_transient(0, '', '');
	TimeParams    = TransientSol.TimeParams;
	loadDirection = TransientSol.DiscreteLoad.direction;

	qProjDispl = zeros(nnMode, TimeParams.numel);
	qProjAccel = zeros(nnMode, TimeParams.numel);

	for inMode = 1:nnMode
		% Get solution for n+1 number of modes.
		[~,  TransientSol] = this_transient(inMode, 'da', '');
		qDispl = TransientSol.ModeDisplacement.q;
		qXDispl = qDispl(nodeList{nodeLabel}.dof(1), :) * loadDirection(1);
		qYDispl = qDispl(nodeList{nodeLabel}.dof(2), :) * loadDirection(2);
		qZDispl = qDispl(nodeList{nodeLabel}.dof(3), :) * loadDirection(3);
		qProjDispl(inMode, :) = qXDispl + qYDispl + qZDispl;
		qAccel= TransientSol.ModeAcceleration.q;
		qXAccel = qAccel(nodeList{nodeLabel}.dof(1), :) * loadDirection(1);
		qYAccel = qAccel(nodeList{nodeLabel}.dof(2), :) * loadDirection(2);
		qZAccel = qAccel(nodeList{nodeLabel}.dof(3), :) * loadDirection(3);
		qProjAccel(inMode, :) = qXAccel + qYAccel + qZAccel;
	end

	% Get the norm of the difference.
	% TODO: See if a mean is more suited (no TimeParams.numel influence).
	qDiffDispl = rms(qProjDispl - qProjDispl(end, :), 2);
	qDiffAccel = rms(qProjAccel - qProjAccel(end, :), 2);

	% Plot the residuals.
	figure("WindowStyle", "docked");
	subplot(1, 2, 1);
	plot_residuals(nModeSet, qDiffDispl, "Mode Displacement");
	subplot(1, 2, 2);
	plot_residuals(nModeSet, qDiffAccel, "Mode Acceleration");

	% TODO: change name
	function plot_residuals(nModeSet, residuals, methodName)
		semilogy(nModeSet, residuals);
		title("Diffs RMS convergence (" + methodName + ")");
		xlabel("Dimension of the modal basis");
		ylabel("Diffs");
		grid;
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

		% TODO: names not evocative.
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
