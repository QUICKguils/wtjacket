function transient_analysis(varargin)
% TRANSIENT_ANALYSIS  Analyze the results of the transient part.
%
% Arguments:
%	focus (1xN char) -- Analysis to focus on (default: 'tscn').
%	  't' -> [T]ransient time response.
%	  's' -> [S]teady-state response.
%	  'c' -> [C]onvergence of the mode superposition methods.
%	  'n' -> [N]ewmark solution.
%	nodeLabels (1xN double) -- Label list of nodes to inspect (default: [18, 22]).

% Set default value for optional inputs.
optargs = {'tscn', [18, 22]};
% Overwrite default value of optional inputs.
optargs(1:numel(varargin)) = varargin;
% Place optional args in memorable variable names.
[focus, nodeLabels] = optargs{:};

% Fetch and overwrite the defaults execution parameters.
RunArg = load_defaults();
RunArg.opts = '';
RunArg.nodeLabels = nodeLabels;
% Fetch the project statement data.
Stm = load_statement();
% Compute the FEM solution.
[~, SdivStruct, AlgSys, FemSol] = modeling(RunArg, Stm);

% Wrapper for the transient response computations.
this_transient = @(nMode, tSet, method) flexible_transient(RunArg, Stm, SdivStruct, AlgSys, FemSol, nMode, tSet, method);

	function varargout = flexible_transient(RunArg, Stm, SdivStruct, AlgSys, FemSol, nMode, tSet, method)
		% FLEXIBLE_TRANSIENT  Allow to use the transient code flexibly.

		% Overwrite the transient execution parameters.
		RunArg.nMode  = nMode;
		RunArg.tSet   = tSet;
		RunArg.method = method;

		% Run the transient code with these execution parameters.
		[AlgSys, TransientSol] = transient(RunArg, Stm, SdivStruct, AlgSys, FemSol);

		% Echo the outputs of transient.
		optrets = {AlgSys, TransientSol};
		varargout(1:nargout) = optrets(1:nargout);
	end

% Analyze the transient code.
if contains(focus, 't'); analyze_transient_response(RunArg, FemSol, this_transient); end
if contains(focus, 's'); analyze_steady_state(RunArg, SdivStruct.nodeList, FemSol, this_transient); end
if contains(focus, 'c'); analyze_convergence(RunArg, this_transient, SdivStruct.nodeList, 1:8); end
if contains(focus, 'n'); analyze_newmark(RunArg, SdivStruct.nodeList, this_transient); end
end

%% 1. Analyze the transient response

function analyze_transient_response(RunArg, FemSol, this_transient)
% ANALYZE_TRANSIENT_RESPONSE
%
% Arguments:
%	RunArg         (struct) -- Transient code execution parameters.
%	FemSol         (struct) -- Solution of the FEM simulation.
%	this_transient (handle) -- Compute the desired transient solution.

% Get the modal superposition parameters.
[~, TransientSol] = this_transient(RunArg.nMode, 0:0.002:10, '');

% Local aliases, for readability.
g   = TransientSol.DiscreteLoad.spatial;
mu  = TransientSol.ModalSup.mu;
x   = FemSol.mode;
eta = TransientSol.ModalSup.eta;
h   = TransientSol.ModalSup.h;
t   = TransientSol.TimeParams.sample;

% 1. Distribution of the modal participation factors

% Amplitude of the modal participation factors.
phiAmpl = abs(x'*g) ./ mu;

figure("WindowStyle", "docked");
stem(phiAmpl);
set(gca,'yscal','log');
grid;
title([ ...
	"Amplitude of the modal participation factors", ...
	"across the first " + num2str(TransientSol.ModalSup.nMode) + " modes"]);
xlabel("Associated mode number");
ylabel("max(\phi(t))");

% 2. Time evolution of the normal coordinates

figure("WindowStyle", "docked");
plot(t, eta);
legend( ...
	["\eta_1", "\eta_2", "\eta_3", "\eta_4", ...
	"\eta_5", "\eta_6", "\eta_7", "\eta_8"], ...
	"NumColumns", 4, ...
	"Location", "south");
title("Time evolution of the normal cooordinates");
xlabel("Time (s)");
ylabel("\eta");
grid;

% 3. Focus on the short period oscillations of mode 5

% Limit the time to the first two seconds.
tCut = t(t<=2);

figure("WindowStyle", "docked");
hold on;
plot(tCut, eta(5, 1:numel(tCut)));
plot(tCut, h(5, 1:numel(tCut)));
legend( ...
	["\eta_5", "h_5"]);
title("Short period oscillations for mode 5");
xlabel("Time (s)");
ylabel("\eta_5, h_5");
grid;
hold off;
end

%% 2. Analyze the steady-state displacement response

function analyze_steady_state(RunArg, nodeList, FemSol, this_transient)
% ANALYZE_STEADY_STATE  Plot the steady-state displacement amplitude.
%
% This function uses the spectral expansion of the admittance matrix to compute
% the amplitude of the displacement.
% It thus assumes that the excitation is harmonic.
%
% Arguments:
%	RunArg         (struct)   -- Transient code execution parameters.
%	nodeList       {1xN Node} -- Cell list of nodes.
%	FemSol         (struct)   -- Solution of the FEM simulation.
%	this_transient (handle)   -- Compute the desired transient solution.

nodeLabel = RunArg.nodeLabels(1);

% Get the modal superposition parameters.
[AlgSys, TransientSol] = this_transient(RunArg.nMode, RunArg.tSet, '');
% Local aliases, for readability.
loadDirection = TransientSol.DiscreteLoad.direction;
g             = TransientSol.DiscreteLoad.spatial;
mu            = TransientSol.ModalSup.mu;

freqHzSet = [ ...
	0    : 0.02   :  0.4,  ...
	0.4  : 0.005  :  0.43, ...
	0.43 : 0.0002 :  0.46, ...
	0.46 : 0.005  :  0.6,  ...
	0.6  : 0.0002 :  0.63, ...
	0.63 : 0.005  :  0.96, ...
	0.96 : 0.0002 :  1,    ...
	1    : 0.005  :  1.3,  ...
	1.3  : 0.1    :  6.8,  ...
	6.8  : 0.02   :  7.6,  ...
	7.6  : 0.2    : 25];
freqRadSet = freqHzSet * 2*pi;

qSteadyAmpl = compute_ss_amplitude(AlgSys, FemSol, mu, g, freqRadSet);

xProj = abs(project_translation(qSteadyAmpl, loadDirection, nodeList, nodeLabel));

figure("WindowStyle", "docked");

semilogy(freqHzSet, xProj);
xlabel("Frequency (Hz)");
ylabel("Displacement (m)");
title( ...
	"Amplitude of the steady harmonic load", ...
	"(node " + num2str(nodeLabel) + ", direction of load)");
grid;

% Find the most interesting frequency for the whales.
[maxDispl, iMaxDisp] = max(xProj);
bestFreqHz = freqHzSet(iMaxDisp);
fprintf( ...
	 ['\nWhales should excite the structure at %.3g Hz\n', ...
	  'in order to maximize the translation at node %u (%.3g m).\n'], ...
	bestFreqHz, nodeLabel, maxDispl);

% Verify that the transient response amplitude converges effectively towards the
% steady-state amplitude that has been found.
%
% Steady-state amplitude, for the load excitation frequency.
xSteady = compute_ss_amplitude(AlgSys, FemSol, mu, g, TransientSol.DiscreteLoad.frequencyRad);
xSteadyProj = abs(project_translation(xSteady, loadDirection, nodeList, nodeLabel));
fprintf( ...
	['\nThe steady-state amplitude of the response at node %u\n', ...
	'under a load frequency of %.3g Hz is %.3g m.\n'], ...
	nodeLabel, TransientSol.DiscreteLoad.frequencyHertz, xSteadyProj);
% Transient response from Newmark method, after a long time.
[~, TransientSol] = this_transient(RunArg.nMode, 0:0.05:100, 'n');
t = TransientSol.TimeParams.sample;
q = project_translation(TransientSol.Newmark.q, loadDirection, nodeList, nodeLabel);
tCut = t((90 <= t) &  (t <= 100));
qCut = q((90 <= t) &  (t <= 100));
figure("WindowStyle", "docked");
plot(tCut, qCut);
xlabel("Time (s)");
ylabel("Displacement (m)");
title( ...
	"Transient response from Newmark method", ...
	"(node " + num2str(nodeLabel) + ", direction of load)");
grid;

	function H = build_H(AlgSys, FemSol, mu, freqRadSet)
		% Discrete set of admittance matrices, for the given frequencies.
		%
		% Arguments:
		%	AlgSys     (struct)         -- Parameters of the discrete algebraic system.
		%	FemSol     (struct)         -- Solution of the FEM simulation.
		%	mu         (nModex1 double) -- Generalized masses.
		%	freqRadSet (1xN double)     -- Harmonic excitation frequencies [rad/s].
		% Return:
		%	H (nDofFreexnDofFreexnMode complex) -- Admittance matrix.

		% Preallocate admittance matrix.
		H = zeros(AlgSys.nDofFree, AlgSys.nDofFree, numel(freqRadSet));
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
		for w = 1:numel(freqRadSet)
			for s = 1:FemSol.nMode
				H(:, :, w) = H(:, :, w) + xTx(:, :, s) ./ (mu(s) * (w0(s)^2 - freqRadSet(w)^2 + 2i*eps(s)*w0(s)*freqRadSet(w)));
			end
		end
	end

	function qSteadyAmpl = compute_ss_amplitude(AlgSys, FemSol, mu, g, freqRadSet)
		% COMPUTE_SS_AMPLITUDE  Compute the steady state response amplitude.
		%
		% This function computes the amplitude of the steady-state displacement,
		% for the given harmonic load frequencies, in [rad/s].

		HSet = build_H(AlgSys, FemSol, mu, freqRadSet);
		qSteadyAmpl = zeros(AlgSys.nDofFree, numel(freqRadSet));
		for w = 1:numel(freqRadSet)
			qSteadyAmpl(:, w) = abs(HSet(:, :, w)) * g;
		end
	end
end

%% 3. Convergence of the mode superposition methods

function analyze_convergence(RunArg, this_transient, nodeList, nModeSet)
% ANALYZE_CONVERGENCE  Convergence of the mode superposition methods.
%
% Arguments:
%	RunArg         (struct)   -- Transient code execution parameters.
%	this_transient (handle)   -- Compute the desired transient solution.
%	nodeList       {1xN Node} -- Cell list of nodes.
%	nModeSet       (1xN int)  -- Set of desired number of first modes.

% NOTE:
% Mode superposition methods require a time step not larger than 0.002s,
% in order to obtain a coherent convergence.
nodeLabel = RunArg.nodeLabels(1);
nnMode    = numel(nModeSet);
tSet      = 0:0.002:10;

% Get the transient solution parameters at disposal.
[~, TransientSol] = this_transient(0, tSet, '');
TimeParams    = TransientSol.TimeParams;
loadDirection = TransientSol.DiscreteLoad.direction;

qProjDispl = zeros(nnMode, TimeParams.numel);
qProjAccel = zeros(nnMode, TimeParams.numel);

for inMode = 1:nnMode
	[~,  TransientSol] = this_transient(nModeSet(inMode), tSet, 'da');
	q = TransientSol.ModeDisplacement.q;
	qProjDispl(inMode, :) = project_translation(q, loadDirection, nodeList, nodeLabel);
	qProjAccel(inMode, :) = project_translation(q, loadDirection, nodeList, nodeLabel);
end

% Reference displacement: Newmark for largest nMode.
[~, TransientSol] = this_transient(nModeSet(end), tSet, 'n');
qProjNewmark = project_translation(TransientSol.Newmark.q, loadDirection, nodeList, nodeLabel);
% Quadratic means of the absolute differences with the reference.
qDiffDispl = rms(qProjDispl - qProjNewmark, 2);
qDiffAccel = rms(qProjAccel - qProjNewmark, 2);

% Plot the absolute differences.
figure("WindowStyle", "docked");
subplot(1, 2, 1);
plot_convergence(nModeSet, qDiffDispl, "mode displacement");
subplot(1, 2, 2);
plot_convergence(nModeSet, qDiffAccel, "mode acceleration");

	function plot_convergence(nModeSet, qDiff, methodName)
		semilogy(nModeSet, qDiff);
		title(["Convergence of RMS of the diffs", "(" + methodName + ")"]);
		xlabel("Dimension of the modal basis");
		ylabel("Diffs");
		grid;
	end
end

%% 4. Analyze the Newmark solution

function analyze_newmark(RunArg, nodeList, this_transient)
% ANALYZE_NEWMARK  Analyze the Newmark solution.
%
% Arguments:
%	RunArg         (struct)   -- Transient code execution parameters.
%	nodeList       {1xN Node} -- Cell list of nodes.
%	this_transient (handle)   -- Compute the desired transient solution.
%
% WARN:
% The fft strongly depend on the chosen time sample.
% It can be advised to use a more refined one than the default.

nodeLabels = RunArg.nodeLabels;

% Get the solution of the transient problem,
% for the Newmark method.
[~, TransientSol] = this_transient(RunArg.nMode, RunArg.tSet, 'n');
% Local aliases, for readability.
timeStep      = TransientSol.TimeParams.steps(1);  % WARN: not robust
loadDirection = TransientSol.DiscreteLoad.direction;
nMode         = TransientSol.ModalSup.nMode;
q             = TransientSol.Newmark.q;

figure("WindowStyle", "docked");

nNode = numel(nodeLabels);
for iNode = 1:nNode
	% Translation along the load direction.
	qProj = project_translation(q, loadDirection, nodeList, nodeLabels(iNode));

	% FFT of the displacement.
	y = fft(qProj);
	% Correct the scaling and shifting.
	fs = 1/timeStep;
	n = length(qProj);
	fshift = (-n/2:n/2-1) * (fs/n);
	yshift = fftshift(y);

	% Plot the FFT.
	subplot(nNode, 1, iNode);
	semilogy(fshift(ceil(end/2:end)), abs(yshift(ceil(end/2:end))));
	title('FFT of the Newmark''s transient response', ...
		['(node: ', num2str(nodeLabels(iNode)), ...
		', order: ', num2str(nMode), ')']);
	xlabel("Frequency (Hz)");
	ylabel("FFT");
	xlim([0, 30]);
	grid;
end
end
