function TimeIntegration = newmark(M, C, K, initialConditions, TimeSet, p)
% NEWMARK  Implement the Newmark's time integration algorithm.
%
% Implement the Newmark's time integration method for linear systems.
%
% Arguments:
%	M, C, K           (nDof x nDof double) -- Matrices of the structural system.
%	initialConditions (nDof x 2 double)    -- Initial conditions.
%	TimeSet           (stuct)              -- 
%	p                 (nDofx1 double)      -- Load applied to the structure [N].
% Return:
%	TimeIntegration (struct) -- Solution of the time integration, with fields:
%	  q       (xxx double) -- Displacements [m].
%	  qd      (xxx double) -- Velocities [m/s].
%	  qdd     (xxx double) -- Accelerations [m/s²].
%	  TimeSet (struct)     -- The TimeSet structure passed in argument.
% This is a straightforward Matlab implementation of the algorithm
% explained in the reference book, section 7.2.1. page 522.

% FIX:
% - Compute for the constrained structure, and not the free one.
% TODO:
% - h is badly implemented.
% PERF:
% - see if the guess are explicit. If so, don't store the in matrices,
% but simple vectors.

% Constant parameters associated with the quadrature scheme.
% See reference book, p. 523.
beta = 0.25;
gamma = 0.5;

% Preallocations.
nDof = size(M, 1);
nCstrDof
q        = zeros(nDof, TimeSet.numel);
qd       = zeros(nDof, TimeSet.numel);
qdd      = zeros(nDof, TimeSet.numel);
q_guess  = zeros(nDof, TimeSet.numel);
qd_guess = zeros(nDof, TimeSet.numel);

% Initial conditions.
q(:, 1)  = initialConditions(:, 1);
qd(:, 1) = initialConditions(:, 2);
%
h = TimeSet.steps;

% Initial accelerations.
qdd(:, 1) = M \ (p(:, 1) - C*qd(:, 1) - K*q(:, 1));

% Iteration matrix.
S = M + gamma*h(1)*C + beta*h(1)^2*K;  % FIX: h(1) not correct

% See reference book, p. 524
for n = 1:TimeSet.numel-1

	% Prediction
	qd_guess(:, n+1) = qd(:, n) + (1-gamma)*h(n)*qdd(:, n);
	q_guess(:, n+1)  = q(:, n) + h(n)*qd(:, n) + (0.5-beta)*h(n)^2*qdd(:, n);

	% Evaluation of accelerations
	qdd(:, n+1) = S \ (p(:, n+1) - C*qd_guess(:, n+1) - K*q_guess(:, n+1));

	% Correction
	qd(:, n+1) = qd_guess(:, n+1) + gamma*h(n)*qdd(:, n+1);
	q(:, n+1)  = q_guess(:, n+1) + beta*h(n)^2*qdd(:, n+1);
end

% Build return data structure.
TimeIntegration.q = q;
TimeIntegration.qdot = qd;
TimeIntegration.qddot = qdd;
TimeIntegration.TimeSet = TimeSet;
TimeIntegration.method = 'Newmark time integration';

end
