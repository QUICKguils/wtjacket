function NewmarkSol = newmark(AlgSys, TimeParams, loadSample)
% NEWMARK  Implement the Newmark's time integration algorithm.
%
% Implement the Newmark's time integration method for linear systems.
%
% Arguments:
%	AlgSys     (struct)                -- Parameters of the discrete algebraic system.
%	TimeParams (struct)                -- Temporal parameters of the problem.
%	loadSample (nDofFreexnTime double) -- Time-discretized load sample [N].
% Return:
%	NewmarkSol (struct)
%	  Solution of the time integration, with fields:
%	    q    (nDofFreexnTime double) -- Displacements [m].
%	    qd   (nDofFreexnTime double) -- Velocities [m/s].
%	    qdd  (nDofFreexnTime double) -- Accelerations [m/s²].
%	    name (1xN char)              -- Name of the method.
%
% This is a straightforward Matlab implementation of the algorithm
% explained in the reference book, section 7.2.1. page 522.

% TODO: h is badly implemented.
% PERF: see if the guess are explicit. If so, don't store them in matrices.

% Constant parameters associated with the quadrature scheme.
beta = 0.25;
gamma = 0.5;

% NOTE:
% As these matrices are only read and not modified, matlab does
% not deep copy them, so these aliases are cheap.
M = AlgSys.M;
C = AlgSys.C;
K = AlgSys.K;

loadSampleCstr = loadSample(~AlgSys.cstrMask, :);

q        = zeros(AlgSys.nDof, TimeParams.numel);
qd       = zeros(AlgSys.nDof, TimeParams.numel);
qdd      = zeros(AlgSys.nDof, TimeParams.numel);
q_guess  = zeros(AlgSys.nDof, TimeParams.numel);
qd_guess = zeros(AlgSys.nDof, TimeParams.numel);

q(:, 1)  = TimeParams.initialConditions(:, 1);
qd(:, 1) = TimeParams.initialConditions(:, 2);

h = TimeParams.steps;

% Initial accelerations.
qdd(:, 1) = M \ (loadSampleCstr(:, 1) - C*qd(:, 1) - K*q(:, 1));

% Iteration matrix.
S = M + gamma*h(1)*C + beta*h(1)^2*K;  % FIX: h(1) not rly correct

for n = 1:TimeParams.numel-1
	% Prediction
	qd_guess(:, n+1) = qd(:, n) + (1-gamma)*h(n)*qdd(:, n);
	q_guess(:, n+1)  = q(:, n)  + h(n)*qd(:, n) + (0.5-beta)*h(n)^2*qdd(:, n);

	% Evaluation of accelerations
	qdd(:, n+1) = S \ (loadSampleCstr(:, n+1) - C*qd_guess(:, n+1) - K*q_guess(:, n+1));

	% Correction
	qd(:, n+1) = qd_guess(:, n+1) + gamma * h(n)   * qdd(:, n+1);
	q(:, n+1)  = q_guess(:, n+1)  + beta  * h(n)^2 * qdd(:, n+1);
end

% Extend the unknown vectors to include the constrained DOFs.
% Anyways, displacements, speeds and accelerations are just null
% for these constrained DOFs.
q_free   = zeros(AlgSys.nDofFree, TimeParams.numel);
qd_free  = zeros(AlgSys.nDofFree, TimeParams.numel);
qdd_free = zeros(AlgSys.nDofFree, TimeParams.numel);
q_free(~AlgSys.cstrMask, :)   = q;
qd_free(~AlgSys.cstrMask, :)  = qd;
qdd_free(~AlgSys.cstrMask, :) = qdd;

NewmarkSol.q     = q_free;
NewmarkSol.qdot  = qd_free;
NewmarkSol.qddot = qdd_free;
NewmarkSol.name  = 'Newmark time integration';
end
