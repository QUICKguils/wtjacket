classdef MainLink < Elem
	% MAINLINK  Represent a main link beam.

	methods
		function ml = MainLink(n1, n2)
			% MAINLINK  Construct an instance of MainLink.

			% Assigned diameter [m] (see project statement).
			diameter = 1;

			ml = ml@Elem(diameter, n1, n2);
		end

		function plotElem(ml, varargin)
			plotElem@Elem(ml, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, varargin{:});
		end
	end
end
