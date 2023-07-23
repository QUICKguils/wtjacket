classdef Elem
	% ELEM Represent a 3D beam element.

	properties
		% Matertial
		rho = 7800;   % Density [kg/m³].
		nu  = 0.3;    % Poisson's ratio [-].
		E   = 210e9;  % Young's modulus [N/m²].
		% Geometry
		w_thick = 0.02;                  % Thickness of the cross-section walls [m].
		d    (1, 1) double {mustBeReal}  % Cross-section diameter [m].
		area (1, 1) double {mustBeReal}  % Cross-section area [m²].
		Iyy  (1, 1) double {mustBeReal}  % Area moment along the local y axis [m^4].
		Izz  (1, 1) double {mustBeReal}  % Area moment along the local z axis [m^4].
		Jx   (1, 1) double {mustBeReal}  % Area moment along the local x axis [m^4].
		% Extremities
		n1 Node  % Starting node.
		n2 Node  % Ending   node.
	end

	% NOTE:
	% Cross-section area and inertia should theoretically be dependant, but
	% these properties need to be inconsistently overridden for the RigidLink
	% class.
	properties (Dependent)
		% Matertial
		mass (1, 1) double {mustBeReal}  % Mass [kg].
		G    (1, 1) double {mustBeReal}  % Shear modulus [N/m²].
		% Geometry
		length (1, 1) double {mustBeReal}  % Length [m].
		dir    (1, 3) double {mustBeReal}  % Directional vector [-].
	end

	methods
		function elem = Elem(d, n1, n2)
			if nargin > 0
				elem.d   = d;
				elem.n1  = n1;
				elem.n2  = n2;

				elem.area = pi * elem.d * elem.w_thick;
				elem.Iyy = pi * ((elem.d+elem.w_thick)^4-(elem.d-elem.w_thick)^4)/64;
				elem.Izz = elem.Iyy;
				elem.Jx = elem.Iyy + elem.Izz;
			end
		end

		function mass = get.mass(elem)
			mass = elem.rho * elem.length * elem.area;
		end

		function G = get.G(elem)
			G = elem.E / (2*(1+elem.nu));
		end

		function length = get.length(elem)
			length = norm(elem.n2.pos - elem.n1.pos);
		end

		function dir = get.dir(elem)
			dir = (elem.n2.pos - elem.n1.pos) / elem.length;
		end

		function plotElem(elem, varargin)
			x = [elem.n1.pos(1), elem.n2.pos(1)];
			y = [elem.n1.pos(2), elem.n2.pos(2)];
			z = [elem.n1.pos(3), elem.n2.pos(3)];
			plot3(x, y, z, varargin{:});
		end
	end
end
