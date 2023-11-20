classdef HarmonicLoad < Load
	% HARMONICLOAD  Represent a harmonic structural load.

	properties (SetAccess = protected)
		amplitude      (1,  1) double {mustBeReal}  % Amplitude [N].
		frequencyHertz (1,  1) double {mustBeReal}  % Frequency [Hz].
		frequencyRad   (1,  1) double {mustBeReal}  % Frequency [rad/s].
	end

	methods
		function hLoad = HarmonicLoad(node, direction, amplitude, frequencyHertz)
			timeEvolution = @(t) amplitude*sin(frequencyHertz*2*pi*t);

			hLoad = hLoad@Load(node, direction, timeEvolution);

			hLoad.amplitude      = amplitude;
			hLoad.frequencyHertz = frequencyHertz;
			hLoad.frequencyRad   = frequencyHertz * 2*pi;
		end

		function DiscreteLoad = set_discrete_load(hLoad, nDofFree, timeSample)
		% SET_DISCRETE_LOAD  Append properties to the superclass method.

			DiscreteLoad = set_discrete_load@Load(hLoad, nDofFree, timeSample);

			% Keep the properties of the HarmonicLoad object.
			DiscreteLoad.amplitude      = hLoad.amplitude;
			DiscreteLoad.frequencyHertz = hLoad.frequencyHertz;
			DiscreteLoad.frequencyRad   = hLoad.frequencyRad;

			% Append the spatial ditribution of the load across the DOFs.
			DiscreteLoad.spatial = zeros(nDofFree, 1);
			DiscreteLoad.spatial(hLoad.node.dof(1)) = hLoad.amplitude * hLoad.direction(1);
			DiscreteLoad.spatial(hLoad.node.dof(2)) = hLoad.amplitude * hLoad.direction(2);
			DiscreteLoad.spatial(hLoad.node.dof(3)) = hLoad.amplitude * hLoad.direction(3);
		end
	end
end
