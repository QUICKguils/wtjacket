classdef AuxLink < Elem
	% AUXLINK  Represent a main link beam.

	methods
		function al = AuxLink(n1, n2)
			% AUXLINK  Construct an instance of AuxLink.

			% Assigned values (see project statement).
			d = 0.6;  % Diameter [m].
			
			al = al@Elem(d, n1, n2);
		end

		function plotElem(al)
			plotElem@Elem(al, Color=[0 0.4470 0.7410], LineWidth=1);
		end
	end
end