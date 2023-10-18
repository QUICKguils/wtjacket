function transient(Constant, ModelingGlobalMatrix, ModelingSolution, SubdivisedStructure, opts)
% TRANSIENT  Transient response due to a collision.

if ModelingSolution.nMode < 8
	warning('wtjacket:NotEnoughModes', ...
		['Not enough modes were computed in modeling part.\n' ...
		'Transient part needs to have at least the first eight modes at its disposal.']);
	return
end

[eps, ~] = set_damping_matrix(Constant.DAMPING_RATIO, ...
	Constant.DAMPING_RATIO, ...
	ModelingSolution.frequencyRad, ...
	ModelingGlobalMatrix.K, ...
	ModelingGlobalMatrix.M);

% TODO: probably put this constant in bare_structure.m or load_constant.m
IMPACT_NODE_LABEL = 18;
TIME_SAMPLE = 0:0.02:10;
load = create_load(Constant, SubdivisedStructure, IMPACT_NODE_LABEL, TIME_SAMPLE);

% TODO: register ICs in Constants.
ModalSuperposition = modal_superposition(ModelingSolution, ModelingGlobalMatrix, eps, load, TIME_SAMPLE, [0, 0]);
q = mode_displacement(ModelingSolution.nMode, ModalSuperposition.nu, ModelingSolution.mode);
plot(TIME_SAMPLE ,q(SubdivisedStructure.nodeList{18}.dof(1), :));
hold on;
plot(TIME_SAMPLE ,q(SubdivisedStructure.nodeList{18}.dof(2), :));
end

function [eps, C] = set_damping_matrix(eps1, eps2, w0, K, M)
% SET_DAMPING_MATRIX  Set the damping matrix, assuming a proportional damping.
%
% Arguments:
%	eps1, eps2 (double)
%	  Damping ratio of the first two modes.
%	frequencies (1 x nMode double)
%	  Natural frequencies of the asociated conservative system, in Hertz.

% Determine the damping ratios for the first eight modes.
% See reference book, p.156.

a = 2             * (w0(1)*eps1 - w0(2)*eps2) / (w0(1)^2 - w0(2)^2);
b = 2*w0(1)*w0(2) * (w0(1)*eps2 - w0(2)*eps1) / (w0(1)^2 - w0(2)^2);

eps = 0.5 * (a*w0 + b./w0);
C   = a*K + b*M;
end


function load = create_load(C, SS, impactNodeLabel, t)
% TODO: make it more general.
% Allow user to enter multiple excited node, force directions and amplitudes.

% Assume an excitation force F = A*sin(w*t).
% TODO: make sure of this chunck of code.
amplitude = 0.5 * C.TAIL_MASS*C.TAIL_SPEED*C.MOMENTUM_TRANSFER / C.IMPACT_DURATION;
frequencyRad = C.LOAD_FREQUENCY * 2*pi;
force  = @(t) amplitude*sin(frequencyRad*t);
forceX = @(t)  force(t) * sind(45);
forceY = @(t) -force(t) * sind(45);

load = zeros(SS.nDof, numel(t));
load(SS.nodeList{impactNodeLabel}.dof(1), :) = forceX(t);
load(SS.nodeList{impactNodeLabel}.dof(2), :) = forceY(t);
end

function ModalSuperposition = modal_superposition(SOL, KM, eps, load, t, initialConditions)
% return the nu and the phi.

% TODO: determine these IC.
% Initial conditions.
A = initialConditions(1);
B = initialConditions(2);

mu  = zeros(SOL.nMode, 1);
phi = zeros(SOL.nMode, length(t));
wd  = zeros(SOL.nMode, 1);
h   = zeros(SOL.nMode, length(t));
nu  = zeros(SOL.nMode, length(t));

for r = 1:SOL.nMode

	mu(r) = SOL.mode(:, r)'*KM.M_free*SOL.mode(:, r);
	% TODO: check mult dims when p(t) is done.
	phi(r, :) = SOL.mode(:, r)' * load / mu(r);

	wd(r) = sqrt(1-eps(r)^2)*SOL.frequencyRad(r);
	h(r, :) = 1/wd(r) .* exp(-eps(r)*wd(r)*t) .* sin(wd(r)*t);
	nu_transient = exp(-eps(r)*wd(r)*t) .* (A*cos(wd(r)*t) + B*sin(wd(r)*t));
	nu_permanent = conv(phi(r, :), h(r, :), 'same');
	% nu_permanent = arrayfun(@(t) integral(@(s) phi{r}(s).*h{r}(t-s), 0, t), tSample);
	nu(r, :) = nu_transient + nu_permanent;
end

ModalSuperposition.phi = phi;
ModalSuperposition.nu  = nu;

end

function q = mode_displacement(k, nu, modes)
% Compute q with the mode displacement method.
q = zeros(size(modes, 1), 1);
for r = 1:k
	q = q + modes(:, r) * nu(r, :);
end
end

function q = mode_acceleration(ModalSuperposition)
% Compute q with the mode acceleration method.
end
