classdef ConcentratedMass
	% CONCENTRATEDMASS  Represent a concentrated mass.

	properties
		% Node location of the concentrated mass.
		node Node
		% Pyhsical parameters.
		mass (1, 1) double {mustBeReal} = 200e3; % Mass [kg].
		Jx   (1, 1) double {mustBeReal} = 24e6;  % Area moment along the X-axis [m^4].
		Iyy  (1, 1) double {mustBeReal} = 24e6;  % Area moment along the Y-axis [m^4].
		Izz  (1, 1) double {mustBeReal} = 24e6;  % Area moment along the Z-axis [m^4].
	end

	methods
		function cm = ConcentratedMass(node, mass, Iyy, Izz, Jx)
			% Mandatory to indicated the location.
			cm.node = node;

			% Other params are optional.
			if nargin > 1
				cm.mass = mass;
				cm.Jx   = Jx;
				cm.Iyy  = Iyy;
				cm.Izz  = Izz;
			end
		end
	end
end