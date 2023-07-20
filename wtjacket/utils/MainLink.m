classdef MainLink < Elem
	% MAINLINK  Represent a main link beam.

	methods
		function ml = MainLink(n1, n2)
			% MAINLINK  Construct an instance of MainLink.

			% Assigned diameter [m] (see project statement).
			d = 1;

			ml = ml@Elem(d, n1, n2);
		end
		
		function plotElem(ml)
			plotElem@Elem(ml, Color=[0.8500 0.3250 0.0980], LineWidth=2);
		end
	end
end