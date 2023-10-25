function transient(C, BS, SS, KM, SOL, opts)
% TRANSIENT  Transient response due to a harmonic excitation.

% 1. Set the time discretization

timeSample = 0:0.001:10;
TimeSet = set_discrete_time(timeSample);

% 2. Set the proportional damping parameters

Damping = set_damping_parameters(C.DAMPING_RATIO, C.DAMPING_RATIO, SOL.frequencyRad, KM);

% 3. Compute the modal superposition

% Choose the load to study.
% NOTE: for the sake of completeness.
% By the way, only one load was assigned in this project.
ThisLoad = BS.loadList{1};
% Create the time-discretized load.
loadSet = ThisLoad.create_load_set(SOL.nDof, TimeSet.sample);

% TODO: this should be passed as transient() argument.
k = SOL.nMode;
ModSup = modal_superposition(k, SOL, KM, Damping.eps, loadSet, TimeSet, C.INITIAL_CONDITIONS);

% 4. Compute the displacements

DisplSetDm = mode_displacement(ModSup.nu, SOL.mode);
DisplSetAm = mode_acceleration(KM, loadSet,  SOL, ModSup);

% 5. Plot the displacements

if contains(opts, 'p')
	lookupNodeLabels = [18, 22];
	plot_displacement(DisplSetDm, TimeSet, lookupNodeLabels, ThisLoad.direction, SS.nodeList, k);
	plot_displacement(DisplSetAm, TimeSet, lookupNodeLabels, ThisLoad.direction, SS.nodeList, k);
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

function ModalSuperposition = modal_superposition(k, SOL, KM, eps, loadSet, TimeSet, initialConditions)
% MODAL_SUPERPOSITION  Compute the modal superposition, from the first vibration modes.

if k > SOL.nMode
	warning("Only the first " + num2str(SOL.nMode) + " modes have been computed.")
	k = SOL.nMode;
end

A = initialConditions(1);
B = initialConditions(2);

mu  = zeros(k, 1);
phi = zeros(k, TimeSet.numel);
wd  = zeros(k, 1);
h   = zeros(k, TimeSet.numel);
nu  = zeros(k, TimeSet.numel);

for r = 1:k
	mu(r)        = SOL.mode(:, r)' * KM.M_free * SOL.mode(:, r);
	phi(r, :)    = SOL.mode(:, r)' * loadSet / mu(r);
	wd(r)        = sqrt(1-eps(r)^2) * SOL.frequencyRad(r);
	h(r, :)      = 1/wd(r) .* exp(-eps(r)*wd(r)*TimeSet.sample) .* sin(wd(r)*TimeSet.sample);
	nuTransient  = exp(-eps(r)*wd(r)*TimeSet.sample) .* (A*cos(wd(r)*TimeSet.sample) + B*sin(wd(r)*TimeSet.sample));
	discreteConv = conv(phi(r, :), h(r, :));
	nuPermanent  = discreteConv(1:length(TimeSet.sample)) .* TimeSet.steps;
	nu(r, :)     = nuTransient + nuPermanent;
end

% TODO: see if rly useful to store h, wd and mu.
ModalSuperposition.phi = phi;
ModalSuperposition.nu  = nu;
ModalSuperposition.mu  = mu;
ModalSuperposition.wd  = wd;
ModalSuperposition.h   = h;
end

%% Compute the displacements

function DisplSet = mode_displacement(nu, mode)
% MODE_DISPLACEMENT  Compute the displacements with the mode displacement method.

DisplSet.q = mode * nu;
DisplSet.method = 'mode displacement';
end

% TODO: this is a non-working draft
function DisplSet = mode_acceleration(KM, loadSet, SOL, MS)
% MODE_ACCELERATION  Compute the displacements with the mode acceleration method.

loadSetCstr = loadSet(~KM.cstrMask, :);
completeStaticResponseCstr = KM.K \ loadSetCstr;
completeStaticResponse = zeros(size(loadSet));
completeStaticResponse(~KM.cstrMask, :) = completeStaticResponseCstr;

partialStaticResponse  = SOL.mode * (MS.phi ./ SOL.frequencyRad.^2);

DisplSet.q = SOL.mode * MS.nu + completeStaticResponse -partialStaticResponse;
DisplSet.method = 'mode acceleration';
end

%% Plot the displacements

function plot_displacement(DisplSet, TimeSet, nodeLabels, dir, nodeList, k)
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
		', order: ', num2str(k), ')']);
	grid;
end
end
