function check_sym(mat, rtol)
% CHECK_SYM Check if a matrix is symmetric, within a certain tolerance.
%   Arguments:
%		mat  (NxN double) -- The matrix to check.
%		rtol (double)     -- The relative tolerance. Default is 1e-6.
%	Throw:
%		(warning) -- Thrown if the symmetry tolerance is not satisfied.

if nargin == 1
	rtol = 1e-6;
end

% WARN: Not the actual relative tolerance.
% The +1 shift is there to avoid divisions by zero.
% This nevertheless gives a decent asymetric check function. Try:
% > fsurf(@(x, y) abs(y-x)/(abs(x)+1))
asym = max(abs(mat'-mat)./(abs(mat)+1), [], "all");
if asym > rtol
	warning('wtjacket:ShouldBeSymmetric', ...
		['Matrix ''' inputname(1) ''' is non-symmetric, but should be. (%g diff.)'], asym);
end

end