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

% TODO: probably put this is bare_structure.m
IMPACT_NODE_LABEL = 18;
p = create_loading(Constant, SubdivisedStructure, IMPACT_NODE_LABEL);

% TODO: register ICs in Constants.
ModalSuperposition = modal_superposition(ModelingSolution, ModelingGlobalMatrix, eps, p, [0, 0]);

q = mode_displacement(8, ModalSuperposition.nu, ModelingSolution.mode);

tSample = 0:0.02:10;
dof = SubdivisedStructure.nodeList{18}.dof(1);  % excitation node, x comp.
qt = q(tSample);
plot(tSample, qt(dof, :));
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


function p = create_loading(C, SS, impactNodeLabel)

% Assume an excitation force F = A*sin(w*t).
A = 0.5 * C.TAIL_MASS*C.TAIL_SPEED*C.MOMENTUM_TRANSFER / C.IMPACT_DURATION;
w = C.LOAD_FREQUENCY * 2*pi;
F  = @(t) A*sin(w*t);
Fx = @(t)  F(t) * sind(45);
Fy = @(t) -F(t) * sind(45);

	% inner function of the closure: isolate parameter `t`.
	function p_at_t = load_at_t(t)
		% In: (double) time instant [s].
		% Out: (1 x nDof double) -- loads for all Dofs at ´t´ [N].

		% PERF: maybe put prealloc outside the inner func
		p_at_t = zeros(SS.nDof, numel(t));
		p_at_t(SS.nodeList{impactNodeLabel}.dof(1), :) = Fx(t);
		p_at_t(SS.nodeList{impactNodeLabel}.dof(2), :) = Fy(t);
	end

	p = @(t) load_at_t(t);
end

function ModalSuperposition = modal_superposition(SOL, KM, eps, p, initialConditions)
% return the nu and the phi.

% TODO: determine these IC.
% Initial conditions.
A = initialConditions(1);
B = initialConditions(2);

mu  = zeros(1, SOL.nMode);
phi = cell(1, SOL.nMode);
wd  = zeros(1, SOL.nMode);
h   = cell(1, SOL.nMode);
nu  = cell(1, SOL.nMode);

% PERF: see if nested function handles are quick enough.
for r = 1:SOL.nMode
	% Generalized masses.
	% TODO: verify if relevant to do that.
	% modes are probably already mass normalized.
	mu(r) = SOL.mode(:, r)'*KM.M_free*SOL.mode(:, r);
	% TODO: check mult dims when p(t) is done.
	phi{r} = @(t) SOL.mode(:, r)' * p(t)/mu(r);

	wd(r) = sqrt(1-eps(r)^2)*SOL.frequencyRad(r);
	h{r} = @(t) 1/wd(r) .* exp(-eps(r)*wd(r)*t) .* sin(wd(r)*t);
	nu_transient = @(t) exp(-eps(r)*wd(r)*t) .* (A*cos(wd(r)*t) + B*sin(wd(r)*t));
	nu_permanent = @(ts) arrayfun(@(t) integral(@(s) phi{r}(s).*h{r}(t-s), 0, t), ts);
	nu{r} = @(t) nu_transient(t) + nu_permanent(t);
end

ModalSuperposition.phi = phi;
ModalSuperposition.nu  = nu;

end

function q = mode_displacement(k, nu, modes)
% Compute q with the mode displacement method.
q = @(t) 0;
for r = 1:k
	q = @(t) q(t) + modes(:, r) * nu{r}(t);
end
end

function q = mode_acceleration(ModalSuperposition)
% Compute q with the mode acceleration method.
end
