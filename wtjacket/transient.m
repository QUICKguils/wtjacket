function [ModalSup, DisplMeth, AccelMeth] = transient(Cst, SdivStruct, FemMat, FemSol, nMode, opts)
% TRANSIENT  Transient response due to a harmonic excitation.
%
% Arguments:
%	Cst        (struct)   -- Constant project quantities.
%	SdivStruct (struct)   -- Subdivised structure.
%	FemMat     (struct)   -- Global structural matrices.
%	FemSol     (struct)   -- Solution of the FEM simulation.
%	nMode      (int)      -- Number of modes involved in the modal superposition.
%	opts       (1xN char) -- Options.
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
% Returns:
%	ModalSup (struct) -- Modal superposition parameters, with fields:
%	  phi (nMode x nTime) -- 
%	  nu  (nMode x nTime) -- 
%	DisplMeth, AccelMeth (struct)
%	  Solution of the mode displacement and mode acceleration method,
%	  with fields:
%	    q      (nDof x nTime) -- Time evolution of the displacements.
%	    method (1 x N char)   -- Name of the method.

% TODO:
% - cleaner to embed TimeSet in the loadSet, DisplMeth and AccelMeth.

% 1. Set the time discretization

timeSample = 0:0.001:10;
TimeSet = set_discrete_time(timeSample);

% 2. Set the proportional damping parameters

Damping = set_damping_parameters(Cst.DAMPING_RATIO, Cst.DAMPING_RATIO, FemSol.frequencyRad, FemMat);

% 3. Compute the modal superposition

% Choose the load to study.
% NOTE: for the sake of completeness.
% By the way, only one load was assigned in this project.
ThisLoad = SdivStruct.loadList{1};
% Create the time-discretized load.
loadSet = ThisLoad.create_load_set(FemSol.nDof, TimeSet.sample);

% TODO: this should be passed as transient() argument.
ModalSup = modal_superposition(nMode, FemSol, FemMat, Damping.eps, loadSet, TimeSet, Cst.INITIAL_CONDITIONS);

% 4. Compute the displacements

DisplMeth = mode_displacement(ModalSup.nu, FemSol.mode, nMode);
AccelMeth = mode_acceleration(FemMat, loadSet,  FemSol, ModalSup, nMode);

% 5. Plot the displacements

if contains(opts, 'p')
	lookupNodeLabels = [18, 22];
	plot_displacement(DisplMeth, TimeSet, lookupNodeLabels, ThisLoad.direction, SdivStruct.nodeList, nMode);
	plot_displacement(AccelMeth, TimeSet, lookupNodeLabels, ThisLoad.direction, SdivStruct.nodeList, nMode);
end

end

%% Set the time discretization

function TimeSet = set_discrete_time(timeSample)
% SET_DISCRETE_TIME  Set the time discretization.

% Ensure the time vector is in the expected shape.
timeSample = reshape(timeSample, 1, []);

% TODO: see if start and end fields are used.
TimeSet.sample = timeSample;
TimeSet.steps  = [diff(timeSample), timeSample(end)-timeSample(end-1)];
TimeSet.start  = timeSample(1);
TimeSet.end    = timeSample(end);
TimeSet.numel  = numel(timeSample);
end

%% Set the proportional damping parameters

function Damping = set_damping_parameters(eps1, eps2, w0, KM)
% SET_DAMPING_MATRIX  Set the damping matrix, assuming a proportional damping.
%
% Arguments:
%	eps1, eps2 (double)
%	  Damping ratio of the first two modes.
%	w0 (1 x nMode double)
%	  Natural frequencies of the asociated conservative system [rad/s].
%	KM (struct)
%	  Global structural matrices.
% Return:
%	Damping (struct)
%	  Damping parameters of the proportional damping model, with fields:
%	    eps (1 x Nmode double) -- Damping ratios of the first modes.
%	    C   (N x N double)     -- Damping matrix.
%
% See reference book, p.156.

a = 2             * (w0(1)*eps1 - w0(2)*eps2) / (w0(1)^2 - w0(2)^2);
b = 2*w0(1)*w0(2) * (w0(1)*eps2 - w0(2)*eps1) / (w0(1)^2 - w0(2)^2);

Damping.eps = 0.5 * (a*w0 + b./w0);
Damping.C   = a*KM.K + b*KM.M;
end

%% Compute the modal superposition

function ModalSup = modal_superposition(nMode, SOL, KM, eps, loadSet, TimeSet, initialConditions)
% MODAL_SUPERPOSITION  Compute the modal superposition, from the first vibration modes.

if nMode > SOL.nMode
	warning("Only the first " + num2str(SOL.nMode) + " modes have been computed.");
	nMode = SOL.nMode;
end

A = initialConditions(1);
B = initialConditions(2);

phi = zeros(nMode, TimeSet.numel);
nu  = zeros(nMode, TimeSet.numel);

for r = 1:nMode
	mu           = SOL.mode(:, r)' * KM.M_free * SOL.mode(:, r);
	phi(r, :)    = SOL.mode(:, r)' * loadSet / mu;
	wd           = sqrt(1-eps(r)^2) * SOL.frequencyRad(r);
	h            = 1/wd .* exp(-eps(r)*wd*TimeSet.sample) .* sin(wd*TimeSet.sample);
	nuTransient  = exp(-eps(r)*wd*TimeSet.sample) .* (A*cos(wd*TimeSet.sample) + B*sin(wd*TimeSet.sample));
	discreteConv = conv(phi(r, :), h);
	nuPermanent  = discreteConv(1:length(TimeSet.sample)) .* TimeSet.steps;
	nu(r, :)     = nuTransient + nuPermanent;
end

ModalSup.phi   = phi;
ModalSup.nu    = nu;
ModalSup.nMode = nMode;
end

%% Compute the displacements

function DisplMeth = mode_displacement(nu, mode, nMode)
% MODE_DISPLACEMENT  Compute the displacements with the mode displacement method.

DisplMeth.q = mode(:, 1:nMode) * nu;
DisplMeth.method = 'mode displacement';
end

% TODO: this is a non-working draft
function AccelMeth = mode_acceleration(FemMat, loadSet, FemSol, ModalSup, nMode)
% MODE_ACCELERATION  Compute the displacements with the mode acceleration method.

loadSetCstr = loadSet(~FemMat.cstrMask, :);
completeStaticResponseCstr = FemMat.K \ loadSetCstr;
completeStaticResponse = zeros(size(loadSet));
completeStaticResponse(~FemMat.cstrMask, :) = completeStaticResponseCstr;

partialStaticResponse  = FemSol.mode(:, 1:nMode) * (ModalSup.phi ./ FemSol.frequencyRad(1:nMode).^2);

AccelMeth.q = FemSol.mode(:, 1:nMode) * ModalSup.nu + completeStaticResponse -partialStaticResponse;
AccelMeth.method = 'mode acceleration';
end

%% Plot the displacements

function plot_displacement(DisplSet, TimeSet, nodeLabels, dir, nodeList, nMode)
% PLOT_DISPLACEMENT  Plot the time evolution of the displacements.

allclose(norm(dir), 1);

figure("WindowStyle", "docked");

nNode = numel(nodeLabels);
for i = 1:nNode
	qX = DisplSet.q(nodeList{nodeLabels(i)}.dof(1), :) * dir(1);
	qY = DisplSet.q(nodeList{nodeLabels(i)}.dof(2), :) * dir(2);
	qZ = DisplSet.q(nodeList{nodeLabels(i)}.dof(3), :) * dir(3);

	qDir = qX + qY + qZ;

	subplot(nNode, 1, i);
	plot(TimeSet.sample, qDir);
	xlabel("Time (s)");
	ylabel("Displacement (dir: [" + num2str(dir, '%.3f  ') + "])");
	title('Transient response', ...
		['(node: ', num2str(nodeLabels(i)), ...
		', method: ', DisplSet.method, ...
		', order: ', num2str(nMode), ')']);
	grid;
end
end
