classdef Elem
	% ELEM Represent a 3D beam element.

	properties (SetAccess = protected)
		% Extremities
		n1 Node  % Starting node.
		n2 Node  % Ending   node.
		% Matertial
		rho = 7800;   % Density [kg/m³].
		E   = 210e9;  % Young's modulus [N/m²].
		nu  = 0.3;    % Poisson's ratio [-].
		% Geometry
		w_thick = 0.02;                     % Thickness of the cross-section walls [m].
		d       (1, 1) double {mustBeReal}  % Cross-section diameter [m].
		area    (1, 1) double {mustBeReal}  % Cross-section area [m²].
		Iyy     (1, 1) double {mustBeReal}  % Area moment along the local y axis [m^4].
		Izz     (1, 1) double {mustBeReal}  % Area moment along the local z axis [m^4].
		Jx      (1, 1) double {mustBeReal}  % Area moment along the local x axis [m^4].
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
		% Elementary local and structural matrices
		K_el (:, :) double {mustBeReal}  % Stiffness [multiple units].
		M_el (:, :) double {mustBeReal}  % Mass      [multiple units].
		K_es (:, :) double {mustBeReal}  % Stiffness [multiple units].
		M_es (:, :) double {mustBeReal}  % Mass      [multiple units].
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

		function K_el = get.K_el(elem)
			% Pre-terms computation for stiffness matrix.
			k = [
				   elem.E*elem.area / elem.length,      elem.G*elem.Jx  / elem.length,   ... % k(1), k(2)
				 2*elem.E*elem.Iyy  / elem.length,    2*elem.E*elem.Izz / elem.length,   ... % k(3), k(4)
				 4*elem.E*elem.Iyy  / elem.length,    4*elem.E*elem.Izz / elem.length,   ... % k(5), k(6)
				 6*elem.E*elem.Iyy  / elem.length^2,  6*elem.E*elem.Izz / elem.length^2, ... % k(7), k(8)
				12*elem.E*elem.Iyy  / elem.length^3, 12*elem.E*elem.Izz / elem.length^3];    % k(9), k(10)

			% Elementary stiffness matrix in local axes.
			K_el = [
				 k(1),      0,     0,     0,     0,     0, -k(1),      0,     0,     0,     0,     0;
				    0,  k(10),     0,     0,     0,  k(8),     0, -k(10),     0,     0,     0,  k(8);
				    0,      0,  k(9),     0, -k(7),     0,     0,      0, -k(9),     0, -k(7),     0;
				    0,      0,     0,  k(2),     0,     0,     0,      0,     0, -k(2),     0,     0;
				    0,      0, -k(7),     0,  k(5),     0,     0,      0,  k(7),     0,  k(3),     0;
				    0,   k(8),     0,     0,     0,  k(6),     0,  -k(8),     0,     0,     0,  k(4);
				-k(1),      0,     0,     0,     0,     0,  k(1),      0,     0,     0,     0,     0;
				    0, -k(10),     0,     0,     0, -k(8),     0,  k(10),     0,     0,     0, -k(8);
				    0,      0, -k(9),     0,  k(7),     0,     0,      0,  k(9),     0,  k(7),     0;
				    0,      0,     0, -k(2),     0,     0,     0,      0,     0,  k(2),     0,     0;
				    0,      0, -k(7),     0,  k(3),     0,     0,      0,  k(7),     0,  k(5),     0;
				    0,   k(8),     0,     0,     0,  k(4),     0,  -k(8),     0,     0,     0, k(6)];
		end

		function M_el = get.M_el(elem)
			% Pre-terms computation for stiffness matrix.
			m = [
				   (elem.d/2)^2 / 3,      (elem.d/2)^2 / 6,   ... % m(1), m(2)
				              1 / 3,                 1 / 6,   ... % m(3), m(4)
				            13 / 35,                9 / 70,   ... % m(5), m(6)
				elem.length*11 / 210, elem.length*13 / 420,   ... % m(7), m(8)
				 elem.length^2 / 105,  elem.length^2 / 140];      % m(9), m(10)

			% Elementary mass matrix in local axes.
			M_el = elem.mass * [
				m(3),     0,     0,    0,      0,      0, m(4),     0,     0,    0,      0,      0;
				   0,  m(5),     0,    0,      0,   m(7),    0,  m(6),     0,    0,      0,  -m(8);
				   0,     0,  m(5),    0,  -m(7),      0,    0,     0,  m(6),    0,   m(8),      0;
				   0,     0,     0, m(1),      0,      0,    0,     0,     0, m(2),      0,      0;
				   0,     0, -m(7),    0,   m(9),      0,    0,     0, -m(8),    0, -m(10),      0;
				   0,  m(7),     0,    0,      0,   m(9),    0,  m(8),     0,    0,      0, -m(10);
				m(4),     0,     0,    0,      0,      0, m(3),     0,     0,    0,      0,      0;
				   0,  m(6),     0,    0,      0,   m(8),    0,  m(5),     0,    0,      0,  -m(7);
				   0,     0,  m(6),    0,  -m(8),      0,    0,     0,  m(5),    0,   m(7),      0;
				   0,     0,     0, m(2),      0,      0,    0,     0,     0, m(1),      0,      0;
				   0,     0,  m(8),    0, -m(10),      0,    0,     0,  m(7),    0,   m(9),      0;
				   0, -m(8),     0,    0,      0, -m(10),    0, -m(7),     0,    0,      0,   m(9)];
		end

		function K_es = get.K_es(elem)
			K_es = elem.local2structural(elem.K_el);
		end

		function M_es = get.M_es(elem)
			M_es = elem.local2structural(elem.M_el);
		end

		function plotElem(elem, varargin)
			x = [elem.n1.pos(1), elem.n2.pos(1)];
			y = [elem.n1.pos(2), elem.n2.pos(2)];
			z = [elem.n1.pos(3), elem.n2.pos(3)];
			plot3(x, y, z, varargin{:});
		end
	end

	methods (Access = protected)
		function smat = local2structural(elem, lmat)
			% LOCAL2STRUCTURAL  From local axes to structural axes.
			%
			% Arguments:
			%	elem (Elem)         -- Structural element.
			%	K_el (12x12 double) -- Stiffness elementary matrix, local axes.
			%	M_el (12x12 double) -- Mass      elementary matrix, local axes.
			% Returns:
			%	K_es (12x12 double) -- Stiffness elementary matrix, structural axes.
			%	M_es (12x12 double) -- Mass      elementary matrix, structural axes.

			% Build the local basis.
			%
			% Normalized x-axis directional vector of the local axes.
			ex = elem.dir';
			% Generate ey and ez by computing the null space of {ex, 0, 0}.
			nullspace = null(ex');
			ey = nullspace(:, 1);
			ez = cross(ex, ey);  % Not nullspace(:,2), to ensure right-handedness.
			lbasis = [ex, ey, ez];

			% Transformation matrix between local and structural axes.
			% The transpose local basis happens to be the rotation operator,
			% as the chosen structural basis is simply the identity matrix.
			T = kron(eye(4), lbasis');

			% Applying the change of basis to K_el and M_el.
			smat = T' * lmat * T;

			check_sym(smat);
		end
	end
end
