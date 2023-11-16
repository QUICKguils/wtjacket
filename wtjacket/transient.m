function varargout = transient(Cst, SdivStruct, AlgSys, FemSol, nMode, opts)
% TRANSIENT  Transient response due to a harmonic excitation.
%
% Arguments:
%	Cst        (struct)   -- Constant project quantities.
%	SdivStruct (struct)   -- Subdivised structure.
%	AlgSys     (struct)   -- Parameters of the discrete algebraic system.
%	FemSol     (struct)   -- Solution of the FEM simulation.
%	nMode      (int)      -- Number of modes used in the modal superposition.
%	opts       (1xN char) -- Options.
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
% Returns:
%	AlgSys   (struct) -- Parameters of the discrete algebraic system, with fields:
%	  K_free    (nDofxnDof double) -- Global siffness matrix, without constraints.
%	  M_free    (nDofxnDof double) -- Global mass matrix, without constraints.
%	  K         (NxN double)       -- Global siffness matrix, with constraints.
%	  M         (NxN double)       -- Global mass matrix, with constraints.
%	  C         (NxN double)       -- Proportional damping matrix, with constraints.
%	  eps       (1xNmode double)   -- Proportional damping ratios.
%	  cstrMask  (1xnDof bool)      -- Index on constrained DOFs.
%	  nDof_free (int)              -- Number of DOFs of the free structure.
%	  nDof      (int)              -- Number of DOFs of the constrained structure.
%	  nCstr     (int)              -- Number of constrained DOFs.
%	TransientSol (struct) -- Solution of the transient problem, with fields:
%	  TimeParams          (struct) -- Temporal parameters of the problem.
%	  DiscreteLoad        (struct) -- time-discretized load.
%	  ModalSup            (struct) -- Modal superposition parameters.
%	  ModeDisplacementSol (struct) -- Solution from the mode displacement method.
%	  ModeAccelerationSol (struct) -- Solution from the mode acceleration method.
%	  NewmarkSol          (struct) -- Solution from the Newmark's time integration.

% TODO:
% - add doc for functions
% - make option to select the desired method.
%   (will be more efficient for convergence study).

% 1. Temporal parameters

% TODO: see if better to pass that in argument of transient().
timeSample = 0:0.01:10;
TimeParams = set_time_parameters(timeSample, Cst.INITIAL_CONDITIONS);

% 2. Proportional damping parameters

eps = [Cst.DAMPING_RATIO, Cst.DAMPING_RATIO];
[AlgSys.C, AlgSys.eps] = set_damping_parameters(eps, FemSol.frequencyRad, AlgSys);

% 3. Time-discretized load

% Choose the load to study.
lookupLoadLabel = 1;
ThisLoad = SdivStruct.loadList{lookupLoadLabel};

% Create the time-discretized load.
DiscreteLoad = ThisLoad.set_discrete_load(AlgSys.nDof_free, TimeParams.sample);

% 4. Compute the modal superposition

ModalSup = modal_superposition(AlgSys, FemSol, TimeParams, DiscreteLoad.sample, nMode);

% 5. Compute the displacements

ModeDisplSol = mode_displacement(FemSol, ModalSup);
ModeAccelSol = mode_acceleration(AlgSys, FemSol, ModalSup, DiscreteLoad.sample);
NewmarkSol   = newmark(AlgSys, TimeParams, DiscreteLoad.sample);

% 6. Plot the displacements

if contains(opts, 'p')
	lookupNodeLabels = [18, 22];
	plot_displacement(ModeDisplSol, TimeParams.sample, lookupNodeLabels, ThisLoad.direction, SdivStruct.nodeList, nMode);
	plot_displacement(ModeAccelSol, TimeParams.sample, lookupNodeLabels, ThisLoad.direction, SdivStruct.nodeList, nMode);
	plot_displacement(NewmarkSol,   TimeParams.sample, lookupNodeLabels, ThisLoad.direction, SdivStruct.nodeList, nMode);
end


% 7. Gather and return the relevant calculated data

TransientSol.TimeParams       = TimeParams;
TransientSol.DiscreteLoad     = DiscreteLoad;
TransientSol.ModalSup         = ModalSup;
TransientSol.ModeDisplacement = ModeDisplSol;
TransientSol.ModeAcceleration = ModeAccelSol;
TransientSol.Newmark          = NewmarkSol;

optrets = {AlgSys, TransientSol};
varargout(1:nargout) = optrets(1:nargout);

end

%% 1. Temporal parameters

function TimeParams = set_time_parameters(timeSample, initialConditions)
% SET_TIME_PARAMETERS  Set the temporal parameters of the problem.

% Ensure the time vector is in the expected shape.
timeSample = reshape(timeSample, 1, []);

% TODO: see if start and end fields are used.
TimeParams.sample = timeSample;
TimeParams.steps  = [diff(timeSample), timeSample(end)-timeSample(end-1)];
TimeParams.start  = timeSample(1);
TimeParams.end    = timeSample(end);
TimeParams.numel  = numel(timeSample);
TimeParams.initialConditions = initialConditions;
end

%% 2. Set the proportional damping parameters

function [C, eps] = set_damping_parameters(eps, w0, AlgSys)
% SET_DAMPING_MATRIX  Set the damping matrix, assuming a proportional damping.
%
% Arguments:
%	eps (1x2 double)
%	  Proportional damping ratios of the first two modes.
%	w0 (1xnMode double)
%	  Natural frequencies of the asociated conservative system [rad/s].
%	AlgSys (struct)
%	  Parameters of the discrete algebraic system.
% Returns:
%	C   (NxN double)     -- Proportional damping matrix.
%	eps (1xNmode double) -- Proportional damping ratios of the first modes.
%
% See reference book, p.156.

a = 2             * (w0(1)*eps(1) - w0(2)*eps(2)) / (w0(1)^2 - w0(2)^2);
b = 2*w0(1)*w0(2) * (w0(1)*eps(2) - w0(2)*eps(1)) / (w0(1)^2 - w0(2)^2);

C   = a*AlgSys.K + b*AlgSys.M;
eps = 0.5 * (a*w0 + b./w0);
end

%% 4. Compute the modal superposition

function ModalSup = modal_superposition(AlgSys, FemSol, TimeParams, loadSample, nMode)
% MODAL_SUPERPOSITION  Compute the modal superposition, from the first vibration modes.
%
% Arguments:
%
% Return:
%	ModalSup (struct) -- Modal superposition parameters, with fields:
%	  mu    (nModex1 double)     -- Generalized masses.
%	  phi   (nModexnTime double) -- Modal participation factors.
%	  nu    (nModexnTime double) -- Modal time functions.
%	  nMode (int)                -- Number of modes used.

if nMode > FemSol.nMode
	warning("Only the first " + num2str(FemSol.nMode) + " modes have been computed.");
	nMode = FemSol.nMode;
end

mu  = zeros(nMode, 1);
phi = zeros(nMode, TimeParams.numel);
eta = zeros(nMode, TimeParams.numel);

% Aliases.
A   = TimeParams.initialConditions(1);
B   = TimeParams.initialConditions(2);
t   = TimeParams.sample;
eps = AlgSys.eps;

for r = 1:nMode
	mu(r)        = FemSol.mode(:, r)' * AlgSys.M_free * FemSol.mode(:, r);
	phi(r, :)    = FemSol.mode(:, r)' * loadSample / mu(r);
	wd           = sqrt(1-eps(r)^2) * FemSol.frequencyRad(r);
	h            = 1/wd .* exp(-eps(r)*wd*t) .* sin(wd*t);
	nuTransient  = exp(-eps(r)*wd*t) .* (A*cos(wd*t) + B*sin(wd*t));
	discreteConv = conv(phi(r, :), h);
	nuPermanent  = discreteConv(1:TimeParams.numel) .* TimeParams.steps;
	eta(r, :)    = nuTransient + nuPermanent;
end

ModalSup.mu    = mu;
ModalSup.phi   = phi;
ModalSup.eta   = eta;
ModalSup.nMode = nMode;
end

%% 5. Compute the displacements

function ModeDisplSol = mode_displacement(FemSol, ModalSup)
% MODE_DISPLACEMENT  Compute the displacements with the mode displacement method.

ModeDisplSol.q    = FemSol.mode(:, 1:ModalSup.nMode) * ModalSup.eta;
ModeDisplSol.name = 'mode displacement';
end

function ModeAccelSol = mode_acceleration(AlgSys, FemSol, ModalSup, loadSample)
% MODE_ACCELERATION  Compute the displacements with the mode acceleration method.

% Mode displacement method.
ModeDisplSol = mode_displacement(FemSol, ModalSup);

% Static response of the complete structure.
loadSampleCstr = loadSample(~AlgSys.cstrMask, :);
completeStaticResponseCstr = AlgSys.K \ loadSampleCstr;
completeStaticResponse = zeros(size(loadSample));
completeStaticResponse(~AlgSys.cstrMask, :) = completeStaticResponseCstr;

% Static response for modes 1 to nMode, already included
% in the mode displacement method.
partialStaticResponse = FemSol.mode(:, 1:ModalSup.nMode) * (ModalSup.phi ./ FemSol.frequencyRad(1:ModalSup.nMode).^2);

ModeAccelSol.q    = ModeDisplSol.q + completeStaticResponse -partialStaticResponse;
ModeAccelSol.name = 'mode acceleration';
end

%% 6. Plot the displacements

function plot_displacement(TransientSol, timeSample, lookupNodeLabels, loadDirection, nodeList, nMode)
% PLOT_DISPLACEMENT  Plot the time evolution of the displacements.

allclose(norm(loadDirection), 1);

figure("WindowStyle", "docked");

nNode = numel(lookupNodeLabels);
for i = 1:nNode
	qX = TransientSol.q(nodeList{lookupNodeLabels(i)}.dof(1), :) * loadDirection(1);
	qY = TransientSol.q(nodeList{lookupNodeLabels(i)}.dof(2), :) * loadDirection(2);
	qZ = TransientSol.q(nodeList{lookupNodeLabels(i)}.dof(3), :) * loadDirection(3);

	qDir = qX + qY + qZ;

	subplot(nNode, 1, i);
	plot(timeSample, qDir);
	xlabel("Time (s)");
	ylabel("Displacement (dir: [" + num2str(loadDirection, '%.3f  ') + "])");
	title('Transient response', ...
		['(node: ', num2str(lookupNodeLabels(i)), ...
		', method: ', TransientSol.name, ...
		', order: ', num2str(nMode), ')']);
	grid;
end
end
