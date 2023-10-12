function allclose(a, b, rtol, atol)
% ALLCLOSE  Returns True if two arrays are element-wise equal within a tolerance.
%
%   Arguments:
%		a, b (NxN double) -- The matrices to check.
%		rtol (double)     -- The relative tolerance. Default is 1e-5.
%		atol (double)     -- The absolute tolerance. Default is 1e-8.
%	Throw:
%		(warning) -- Thrown if the symmetry tolerance is not satisfied.
%
%	NOTE:
%	This function is inspired from the numpy `allclose()` function. See:
%	https://numpy.org/doc/stable/reference/generated/numpy.allclose.html

if nargin == 2
	rtol = 1e-5;
	atol = 1e-8;
end

if abs(a - b) > (atol + rtol * abs(b))
	% Identify the maximum absolute difference and its location.
	[diff, ind] = max(abs(a - b), [], 'all');
	[row, col] = ind2sub(size(a), ind);
	% Throw a warning.
	warning( ...
		'wtjacket:ShouldBeSymmetric', ...
		['Matrix ''' inputname(1) ''' and ''' inputname(2) ''' are not close, but should be.\n' ...
			'(diff: %g at (%g, %g).)'], ...
		diff, row, col);
end
end