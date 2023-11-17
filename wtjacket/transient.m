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

%% 4. Compute the modal superposition

function ModalSup = modal_superposition(AlgSys, FemSol, TimeParams, loadSample, nMode)
% MODAL_SUPERPOSITION  Compute the modal superposition, from the first vibration modes.
%
% Arguments:
%
% Return:
%	ModalSup (struct) -- Modal superposition parameters, with fields:
%	  phi   (nModexnTime double) -- Modal participation factors.
%	  nu    (nModexnTime double) -- Modal time functions.
%	  nMode (int)                -- Number of modes used.

if nMode > FemSol.nMode
	warning("Only the first " + num2str(FemSol.nMode) + " modes have been computed.");
	nMode = FemSol.nMode;
end

phi = zeros(nMode, TimeParams.numel);
eta = zeros(nMode, TimeParams.numel);

% Aliases.
A   = TimeParams.initialConditions(1);
B   = TimeParams.initialConditions(2);
t   = TimeParams.sample;
eps = AlgSys.eps;

for r = 1:nMode
	mu           = FemSol.mode(:, r)' * AlgSys.M_free * FemSol.mode(:, r);
	phi(r, :)    = FemSol.mode(:, r)' * loadSample / mu;
	wd           = sqrt(1-eps(r)^2) * FemSol.frequencyRad(r);
	h            = 1/wd .* exp(-eps(r)*wd*t) .* sin(wd*t);
	nuTransient  = exp(-eps(r)*wd*t) .* (A*cos(wd*t) + B*sin(wd*t));
	discreteConv = conv(phi(r, :), h);
	nuPermanent  = discreteConv(1:TimeParams.numel) .* TimeParams.steps;
	eta(r, :)    = nuTransient + nuPermanent;
end

ModalSup.phi        = phi;
ModalSup.eta        = eta;
ModalSup.nMode      = nMode;
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