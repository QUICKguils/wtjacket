classdef Load
	% LOAD  Represent a harmonic structural load.

	properties (SetAccess = protected)
		node           Node                         % Node at which the load is applied.
		direction      (1,  3) double {mustBeReal}  % Direction (structural axes).
		amplitude      (1,  1) double {mustBeReal}  % Amplitude [N].
		frequencyHertz (1,  1) double {mustBeReal}  % Frequency [Hz].
		frequencyRad   (1,  1) double {mustBeReal}  % Frequency [rad/s].
	end

	methods
		function load = Load(node, direction, amplitude, frequencyHertz)
			if nargin > 0
				load.node           = node;
				load.direction      = direction;
				load.amplitude      = amplitude;
				load.frequencyHertz = frequencyHertz;
				load.frequencyRad   = frequencyHertz * 2*pi;
			end
		end

		function loadSet = create_load_set(load, nDof, timeSample)
		% CREATE_LOAD_SET  Create a set of time-discretized loads.
		%
		% Arguments:
		%	nDof       (double)       -- Number of structural DOFs.
		%	timeSample (1 x N double) -- Time sample [s].

			harmonicLoad = @(t) load.amplitude*sin(load.frequencyRad*t);
			harmonicLoadX= @(t)  harmonicLoad(t) * load.direction(1);
			harmonicLoadY= @(t)  harmonicLoad(t) * load.direction(2);
			harmonicLoadZ= @(t)  harmonicLoad(t) * load.direction(3);

			loadSet = zeros(nDof, numel(timeSample));
			loadSet(load.node.dof(1), :) = harmonicLoadX(timeSample);
			loadSet(load.node.dof(2), :) = harmonicLoadY(timeSample);
			loadSet(load.node.dof(3), :) = harmonicLoadZ(timeSample);
		end
	end
end
