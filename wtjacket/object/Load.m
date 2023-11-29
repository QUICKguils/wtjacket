classdef Load
	% LOAD  Represent a structural load.

	properties (SetAccess = protected)
		node          Node                        % Node at which the load is applied.
		direction     (1, 3) double {mustBeReal}  % Direction (structural axes) [m].
		timeEvolution function_handle             % Time evolution law [s -> N].
	end

	methods
		function load = Load(node, direction, timeEvolution)
			if nargin > 0
				load.node          = node;
				load.direction     = direction;
				load.timeEvolution = timeEvolution;

			end
		end

		function DiscreteLoad = set_discrete_load(load, nDofFree, timeSample)
			% SET_DISCRETE_LOAD  Create a sample of time-discretized load.
			%
			% Arguments:
			%	nDofFree  (int)             -- Number of DOFs of the free structure.
			%	timeSample (1xnTime double) -- Time sample [s].
			% Return:
			%	DiscreteLoad (struct) -- time-discretized load, with fields:
			%	  node          (Node)              -- Node at which the load is applied.
			%	  direction     (1x3 double)        -- Direction (structural axes) [m].
			%	  timeEvolution (f_handle)          -- Time evolution law [s -> N].
			%	  sample        (nDofxnTime double) -- Load sample [N].

			% Keep the properties of the Load object.
			DiscreteLoad.node          = load.node;
			DiscreteLoad.direction     = load.direction;
			DiscreteLoad.timeEvolution = load.timeEvolution;

			load_X = @(t) load.timeEvolution(t) * load.direction(1);
			load_Y = @(t) load.timeEvolution(t) * load.direction(2);
			load_Z = @(t) load.timeEvolution(t) * load.direction(3);

			DiscreteLoad.sample = zeros(nDofFree, numel(timeSample));
			DiscreteLoad.sample(load.node.dof(1), :) = load_X(timeSample);
			DiscreteLoad.sample(load.node.dof(2), :) = load_Y(timeSample);
			DiscreteLoad.sample(load.node.dof(3), :) = load_Z(timeSample);
		end
	end
end
