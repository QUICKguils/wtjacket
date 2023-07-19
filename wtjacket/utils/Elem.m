classdef Elem
	% ELEM Represent a 3D beam element.

	properties
		% Matertial properties.
		rho = 7800;                     % Density [kg/m³].
		nu  = 0.3;                      % Poisson's ratio [-].
		E   = 210e9;                    % Young's modulus [N/m²].
		G   (1, 1) double {mustBeReal}  % Shear modulus [N/m²].
		% Geometric properties.
		d   (1, 1) double {mustBeReal}  % Cross-section diameter [m].
		A   (1, 1) double {mustBeReal}  % Cross-section area [m²].
		Iyy (1, 1) double {mustBeReal}  % Area moment along the local y axis [m^4].
		Izz (1, 1) double {mustBeReal}  % Area moment along the local z axis [m^4].
		Jx  (1, 1) double {mustBeReal}  % Area moment along the local x axis [m^4].
		% Extremities.
		n1 Node  % Starting node.
		n2 Node  % Ending   node.
	end

	methods
		function elem = Elem(d, n1, n2, varargin)
			% ELEM  Construct an instance of Elem.
			if nargin > 0
				elem.d   = d;
				elem.n1  = n1;
				elem.n2  = n2;

				elem.G = elem.E / (2*(1+elem.nu));
				elem.A   = 2*pi * elem.d^2/4;
				elem.Iyy = pi * (6*elem.d^3+8*elem.d)/64;  % tubular section
				elem.Izz = elem.Iyy;
				elem.Jx  = elem.Iyy + elem.Izz;
			end
		end

		function l = length(elem)
			% LENGTH  Compute the length of the elem.
			l = norm(elem.n2.pos - elem.n1.pos);
		end

		function plotElem(elem, varargin)
			x = [elem.n1.pos(1), elem.n2.pos(1)];
			y = [elem.n1.pos(2), elem.n2.pos(2)];
			z = [elem.n1.pos(3), elem.n2.pos(3)];
			plot3(x, y, z, varargin{:});
		end
	end
end
