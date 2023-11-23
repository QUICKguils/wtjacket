function varargout = transient(Cst, SdivStruct, AlgSys, FemSol, nMode, method, opts)
% TRANSIENT  Transient response due to a harmonic excitation.
%
% Arguments:
%	Cst        (struct)   -- Constant project quantities.
%	SdivStruct (struct)   -- Subdivised structure.
%	AlgSys     (struct)   -- Parameters of the discrete algebraic system.
%	FemSol     (struct)   -- Solution of the FEM simulation.
%	nMode      (int)      -- Number of modes used in the modal superposition.
%	method     (1xN char) -- Methods used to compute the transient response.
%	  'd' -> Mode [D]isplacement method.
%	  'a' -> Mode [A]cceleration method.
%	  'n' -> [N]ewmark (time integration).
%	opts       (1xN char) -- Options.
%	  'p' -> Enable [P]lots creation.
% Returns:
%	AlgSys (struct) -- Parameters of the discrete algebraic system, with fields:
%	  K_free   (nDofFreexnDofFree double) -- Global siffness matrix, without constraints.
%	  M_free   (nDofFreexnDofFree double) -- Global mass matrix, without constraints.
%	  K        (nDofxnDof double)         -- Global siffness matrix, with constraints.
%	  M        (nDofxnDof double)         -- Global mass matrix, with constraints.
%	  C        (nDofxnDof double)         -- Proportional damping matrix, with constraints.
%	  eps      (1xNmode double)           -- Proportional damping ratios.
%	  cstrMask (1xnDofFree bool)          -- Index on constrained DOFs.
%	  nDofFree (int)                      -- Number of DOFs of the free structure.
%	  nDof     (int)                      -- Number of DOFs of the constrained structure.
%	  nCstr    (int)                      -- Number of constrained DOFs.
%	TransientSol (struct) -- Solutions of the transient problem, with fields:
%	  TimeParams       (struct) -- Temporal parameters of the problem.
%	  DiscreteLoad     (struct) -- time-discretized load.
%	  ModalSup         (struct) -- Modal superposition parameters.
%	  ModeDisplacement (struct) -- Solution from the mode displacement method.
%	  ModeAcceleration (struct) -- Solution from the mode acceleration method.
%	  Newmark          (struct) -- Solution from the Newmark's time integration.

% TODO:
% - carry conditional logic on method in this main function.

% 1. Temporal parameters

% TODO:
% Allow time sample to be defined for each methods.
% Newmark does not really require such a small time step
% to give satifying results.
timeSample = 0:0.01:10;
TransientSol.TimeParams = set_time_parameters(timeSample, Cst.INITIAL_CONDITIONS);

% 2. Proportional damping parameters

eps = [Cst.DAMPING_RATIO, Cst.DAMPING_RATIO];
[AlgSys.C, AlgSys.eps] = set_damping_parameters(eps, FemSol.frequencyRad, AlgSys);

% 3. Time-discretized load

% Choose the load to study.
inspectLoadLabel = 1;
ThisLoad = SdivStruct.loadList{inspectLoadLabel};

% Create the time-discretized load.
TransientSol.DiscreteLoad = ThisLoad.set_discrete_load(AlgSys.nDofFree, TransientSol.TimeParams.sample);

% 4. Compute the modal superposition

TransientSol.ModalSup = modal_superposition(AlgSys, FemSol, TransientSol, nMode);

% 5. Compute the displacements

TransientSol = compute_displacement(AlgSys, FemSol, TransientSol, method);

% 6. Plot the displacements

if contains(opts, 'p')
	inspectNodeLabels = [18, 22];
	plot_displacement(TransientSol, inspectNodeLabels, SdivStruct.nodeList, method);
end

% 7. Return the relevant calculated data

optrets = {AlgSys, TransientSol};
varargout(1:nargout) = optrets(1:nargout);

end

%% 1. Temporal parameters

function TimeParams = set_time_parameters(timeSample, initialConditions)
% SET_TIME_PARAMETERS  Set the temporal parameters of the problem.
%
% Arguments:
%	timeSample        (1xN double) -- Time sample [s].
%	initialConditions (1x2 double) -- System's initial conditions.

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

function ModalSup = modal_superposition(AlgSys, FemSol, TransientSol, nMode)
% MODAL_SUPERPOSITION  Compute the modal superposition, from the first vibration modes.
%
% Arguments:
%	AlgSys       (struct) -- Parameters of the discrete algebraic system.
%	FemSol       (struct) -- Solution of the FEM simulation.
%	TransientSol (struct) -- Parameters of the transient problem.
%	nMode        (int)    -- Number of modes to use.
% Return:
%	ModalSup (struct) -- Modal superposition parameters, with fields:
%	  mu    (nModex1 double)     -- Generalized masses.
%	  wd    (nModex1 double)     -- Damped natural frequencies [rad/s].
%	  phi   (nModexnTime double) -- Modal participation factors.
%	  nu    (nModexnTime double) -- Modal time functions.
%	  nMode (int)                -- Number of modes used.

if nMode > FemSol.nMode
	warning("Only the first " + num2str(FemSol.nMode) + " modes have been computed.");
	nMode = FemSol.nMode;
end

TimeParams = TransientSol.TimeParams;
loadSample = TransientSol.DiscreteLoad.sample;
A          = TimeParams.initialConditions(1);
B          = TimeParams.initialConditions(2);
t          = TimeParams.sample;
eps        = AlgSys.eps;

mu  = zeros(nMode, 1);
wd  = zeros(nMode, 1);
phi = zeros(nMode, TimeParams.numel);
eta = zeros(nMode, TimeParams.numel);

for r = 1:nMode
	mu(r)        = FemSol.mode(:, r)' * AlgSys.M_free * FemSol.mode(:, r);
	phi(r, :)    = FemSol.mode(:, r)' * loadSample / mu(r);
	wd(r)        = sqrt(1-eps(r)^2) * FemSol.frequencyRad(r);
	h            = 1/wd(r) .* exp(-eps(r)*wd(r)*t) .* sin(wd(r)*t);
	nuTransient  = exp(-eps(r)*wd(r)*t) .* (A*cos(wd(r)*t) + B*sin(wd(r)*t));
	discreteConv = conv(phi(r, :), h);
	nuPermanent  = discreteConv(1:TimeParams.numel) .* TimeParams.steps;
	eta(r, :)    = nuTransient + nuPermanent;
end

ModalSup.mu    = mu;
ModalSup.wd    = wd;
ModalSup.phi   = phi;
ModalSup.eta   = eta;
ModalSup.nMode = nMode;
end

%% 5. Compute the displacements

function ModeDisplacement = mode_displacement(FemSol, ModalSup)
% MODE_DISPLACEMENT  Compute the displacements with the mode displacement method.

ModeDisplacement.q    = FemSol.mode(:, 1:ModalSup.nMode) * ModalSup.eta;
ModeDisplacement.name = 'mode displacement';
end

function ModeAcceleration = mode_acceleration(AlgSys, FemSol, TransientSol)
% MODE_ACCELERATION  Compute the displacements with the mode acceleration method.

ModalSup   = TransientSol.ModalSup;
loadSample = TransientSol.DiscreteLoad.sample;

% Mode displacement method.
ModeDisplMethod = mode_displacement(FemSol, ModalSup);

% Static response of the complete structure.
loadSampleCstr = loadSample(~AlgSys.cstrMask, :);
completeStaticResponseCstr = AlgSys.K \ loadSampleCstr;
completeStaticResponse = zeros(size(loadSample));
completeStaticResponse(~AlgSys.cstrMask, :) = completeStaticResponseCstr;

% Static response for modes 1 to nMode, already included
% in the mode displacement method.
partialStaticResponse = FemSol.mode(:, 1:ModalSup.nMode) * (ModalSup.phi ./ FemSol.frequencyRad(1:ModalSup.nMode).^2);

ModeAcceleration.q    = ModeDisplMethod.q + completeStaticResponse -partialStaticResponse;
ModeAcceleration.name = 'mode acceleration';
end

function TransientSol = compute_displacement(AlgSys, FemSol, TransientSol, method)
% COMPUTE_DISPLACEMENT  Select the methods for which transient solutions are desired.

ModalSup     = TransientSol.ModalSup;
TimeParams   = TransientSol.TimeParams;
DiscreteLoad = TransientSol.DiscreteLoad;

if contains(method, 'd'); TransientSol.ModeDisplacement = mode_displacement(FemSol, ModalSup); end
if contains(method, 'a'); TransientSol.ModeAcceleration = mode_acceleration(AlgSys, FemSol, TransientSol); end
if contains(method, 'n'); TransientSol.Newmark          = newmark(AlgSys, TimeParams, DiscreteLoad.sample); end
end

%% 6. Plot the displacements

function plot_one_displacement(TransientSol, Method, inspectNodeLabels, nodeList)
% PLOT_ONE_DISPLACEMENT  Plot the displacements from one given method.
%
% Arguments:
%	TransientSol      (struct)   -- Solutions of the transient problem.
%	Method            (struct)   -- Solution from one transient method.
%	inspectNodeLabels (1xN int)  -- Labels list of node to inspect.
%	nodeList          {1xN Node} -- Cell list of nodes.

timeSample    = TransientSol.TimeParams.sample;
loadDirection = TransientSol.DiscreteLoad.direction;
nMode         = TransientSol.ModalSup.nMode;

allclose(norm(loadDirection), 1);

figure("WindowStyle", "docked");

nNode = numel(inspectNodeLabels);
for iNode = 1:nNode
	qX = Method.q(nodeList{inspectNodeLabels(iNode)}.dof(1), :) * loadDirection(1);
	qY = Method.q(nodeList{inspectNodeLabels(iNode)}.dof(2), :) * loadDirection(2);
	qZ = Method.q(nodeList{inspectNodeLabels(iNode)}.dof(3), :) * loadDirection(3);

	qDir = qX + qY + qZ;

	subplot(nNode, 1, iNode);
	plot(timeSample, qDir);
	xlabel("Time (s)");
	ylabel("Displacement (dir: [" + num2str(loadDirection, '%.3f  ') + "])");
	title('Transient response', ...
		['(node: ', num2str(inspectNodeLabels(iNode)), ...
		', method: ', Method.name, ...
		', order: ', num2str(nMode), ')']);
	grid;
end
end

function plot_displacement(TransientSol, inspectNodeLabels, nodeList, method)
% PLOT_DISPLACEMENT  Select the methods for which displacement plots are desired.
%
% Arguments:
%	TransientSol      (struct)   -- Solutions of the transient problem.
%	inspectNodeLabels (1xN int)  -- Labels list of node to inspect.
%	nodeList          {1xN Node} -- Cell list of nodes.
%	method            (1xN char) -- Methods used to compute the transient response.
%	  'd' -> Mode displacement method.
%	  'a' -> Mode acceleration method.
%	  'n' -> Newmark (time integration).

plot_for_method = @(Method) plot_one_displacement(TransientSol, Method, inspectNodeLabels, nodeList);

if contains(method, 'd'); plot_for_method(TransientSol.ModeDisplacement); end
if contains(method, 'a'); plot_for_method(TransientSol.ModeAcceleration); end
if contains(method, 'n'); plot_for_method(TransientSol.Newmark); end
end
