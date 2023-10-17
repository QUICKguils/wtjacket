function transient(C, KM, SOL, opts)
% TRANSIENT  Transient response due to a collision.

if SOL.nMode < 8
	warning('wtjacket:NotEnoughModes', ...
		['Not enough modes were computed in modeling part.\n' ...
		'Transient part needs to have at least the first eight modes at its disposal.']);
	return
end

[eps, dampingMatrix] = set_damping_matrix(C.DAMPING_RATIO, C.DAMPING_RATIO, SOL.frequencies, KM.K, KM.M);
disp(eps);

end

function [eps, dampingMatrix] = set_damping_matrix(eps1, eps2, frequencies, K, M)
% SET_DAMPING_MATRIX  Set the damping matrix, assuming a proportional damping.
%
% Arguments:
%	eps1, eps2 (double)
%	  Damping ratio of the first two modes.
%	frequencies (1 x nMode double)
%	  Natural frequencies of the asociated conservative system, in Hertz.

% Determine the damping ratios for the first eight modes.
% See reference book, p.156.
w = frequencies * 2*pi;

a = 2           * (w(1)*eps1 - w(2)*eps2) / (w(2)^1 - w(1)^2);
b = 2*w(1)*w(2) * (w(2)*eps1 - w(1)*eps2) / (w(2)^2 - w(2)^1);

eps = 0.5 * (a*w + b./w);
dampingMatrix = a*K + b*M;
end