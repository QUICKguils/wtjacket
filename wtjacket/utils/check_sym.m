function check_sym(mat, atol)
% CHECK_SYM Check if a matrix is symmetric, within a certain tolerance.
%   Arguments:
%		mat  (NxN double) -- The matrix to check.
%		atol (double)     -- The absolute tolerance. Default is 1e-4.
%	Throw:
%		(error) -- Thrown if the symmetry tolerance is not satisfied.

if nargin == 1
	atol = 1e-4;
end

asym = max(mat-mat', [], "all");
assert(asym < atol, ...
	'Matrix is non-symmetric, but should be. (%g diff.)', asym);
end