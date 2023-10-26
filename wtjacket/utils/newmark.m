function NewmarkSol = newmark(M, C, K, TimeParams, loadSample)
% NEWMARK  Implement the Newmark's time integration algorithm.
%
% Implement the Newmark's time integration method for linear systems.
%
% Arguments:
%	M, C, K    (nDofxnDof double)  -- Matrices of the structural system.
%	TimeParams (stuct)             -- Temporal parameters of the problem.
%	loadSample (nDofxnTime double) -- time-discretized load sample [N].
% Return:
%	NewmarkSol (struct) -- Solution of the time integration, with fields:
%	  q          (nDofxnTime double) -- Displacements [m].
%	  qd         (nDofxnTime double) -- Velocities [m/s].
%	  qdd        (nDofxnTime double) -- Accelerations [m/s²].
%	  timeSample (1xnTime)           -- Time sample used [s].
%	  loadSample (nDofxnTime double) -- time-discretized load sample [N].
%
% This is a straightforward Matlab implementation of the algorithm
% explained in the reference book, section 7.2.1. page 522.

% FIX: Compute for the constrained structure, and not the free one.
% TODO: h is badly implemented.
% PERF: see if the guess are explicit. If so, don't store them in matrices, but simple vectors.

% Constant parameters associated with the quadrature scheme.
beta = 0.25;
gamma = 0.5;

nDof     = size(M, 1);
q        = zeros(nDof, TimeParams.numel);
qd       = zeros(nDof, TimeParams.numel);
qdd      = zeros(nDof, TimeParams.numel);
q_guess  = zeros(nDof, TimeParams.numel);
qd_guess = zeros(nDof, TimeParams.numel);

q(:, 1)  = TimeParams.initialConditions(:, 1);
qd(:, 1) = TimeParams.initialConditions(:, 2);

h = TimeParams.steps;

% Initial accelerations.
qdd(:, 1) = M \ (loadSample(:, 1) - C*qd(:, 1) - K*q(:, 1));

% Iteration matrix.
S = M + gamma*h(1)*C + beta*h(1)^2*K;  % FIX: h(1) not rly correct

for n = 1:TimeParams.numel-1
	% Prediction
	qd_guess(:, n+1) = qd(:, n) + (1-gamma)*h(n)*qdd(:, n);
	q_guess(:, n+1)  = q(:, n) + h(n)*qd(:, n) + (0.5-beta)*h(n)^2*qdd(:, n);

	% Evaluation of accelerations
	qdd(:, n+1) = S \ (loadSample(:, n+1) - C*qd_guess(:, n+1) - K*q_guess(:, n+1));

	% Correction
	qd(:, n+1) = qd_guess(:, n+1) + gamma*h(n)*qdd(:, n+1);
	q(:, n+1)  = q_guess(:, n+1) + beta*h(n)^2*qdd(:, n+1);
end

NewmarkSol.q          = q;
NewmarkSol.qdot       = qd;
NewmarkSol.qddot      = qdd;
NewmarkSol.timeSample = TimeParams.sample;
NewmarkSol.loadSample = loadSample;
NewmarkSol.name       = 'Newmark time integration';
end
