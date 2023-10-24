function TimeIntegration = newmark(M, C, K, initialConditions, TimeSet, p)
% NEWMARK  Implement the Newmark's time integration algorithm.

% Constant parameters associated with the quadrature scheme.
% See reference book, p. 523.
beta = 0.25;
gamma = 0.5;

% Preallocations.
nDof = size(M, 1);
q          = zeros(nDof, TimeSet.numel);
qdot       = zeros(nDof, TimeSet.numel);
qddot      = zeros(nDof, TimeSet.numel);
q_guess    = zeros(nDof, TimeSet.numel);
qdot_guess = zeros(nDof, TimeSet.numel);

% Initial conditions.
q(:, 1)    = initialConditions(:, 1);
qdot(:, 0) = initialConditions(:, 2);
%
h = TimeSet.step;

% Initial accelerations.
qddot(:, 1) = M \ (p(:, 1) - C*qdot_0 - K*q_0);

% Iteration matrix.
S = M + gamma*h*C + beta*h^2*K;

% See reference book, p. 524
for t = 1:TimeSet.numel-1

	% Prediction
	qdot_guess(:, t+1) = qdot(:, t) + (1-gamma)*h(t)*qddot(:, t);
	q_guess(:, t+1)    = q_n + h*qdot(:, t) + (0.5-beta)*h(t)^2*qddot(:, t);

	% Evaluation of accelerations
	qddot(:, t+1) = S \ (p(:, t+1) - C*qdot_guess(:, t+1) - K*q_guess(:, t+1));

	% Correction
	qdot(:, t+1) = qdot_guess_n1 + gamma*h(t)*qddot(:, t+1);
	q(:, t+1)    = q_guess(:, t+1) + beta*h^2*qddot(:, t+1);
end

% Build return data structure.
TimeIntegration.q = q;
TimeIntegration.qdot = qdot;
TimeIntegration.qddot = qddot;
TimeIntegration.TimeSet = TimeSet;

end