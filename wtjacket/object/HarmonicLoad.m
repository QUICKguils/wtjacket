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

		function DiscreteLoad = set_discrete_load(load, nDof, timeSample)
		% SET_DISCRETE_LOAD  Append properties to the superclass method.

			DiscreteLoad = set_discrete_load@Load(load, nDof, timeSample);

			% Keep the properties of the HarmonicLoad object.
			DiscreteLoad.amplitude      = load.amplitude;
			DiscreteLoad.frequencyHertz = load.frequencyHertz;
			DiscreteLoad.frequencyRad   = load.frequencyRad;

			% Append the spatial ditribution of the load across the DOFs.
			DiscreteLoad.spatial = zeros(nDof, 1);
			DiscreteLoad.spatial(load.node.dof(1)) = load.amplitude * load.direction(1);
			DiscreteLoad.spatial(load.node.dof(2)) = load.amplitude * load.direction(2);
			DiscreteLoad.spatial(load.node.dof(3)) = load.amplitude * load.direction(3);
		end
	end
end
