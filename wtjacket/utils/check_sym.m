function check_sym(mat, atol)
% CHECK_SYM Check if a matrix is symmetric, within a certain tolerance.
%   Arguments:
%		mat  (NxN double) -- The matrix to check.
%		atol (double)     -- The absolute tolerance. Default is 1e-4.
%	Throw:
%		(warning) -- Thrown if the symmetry tolerance is not satisfied.

if nargin == 1
	atol = 1e-4;
end

asym = max(mat-mat', [], "all");
if asym > atol
	warning(['Matrix ''' inputname(1) ''' is non-symmetric, but should be. (%g diff.)'], asym);
end