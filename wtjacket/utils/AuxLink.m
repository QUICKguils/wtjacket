classdef AuxLink < Elem
	% AUXLINK  Represent a main link beam.

	methods
		function al = AuxLink(n1, n2)
			% AUXLINK  Construct an instance of AuxLink.

			% Assigned diameter [m] (see project statement).
			d = 0.6;

			al = al@Elem(d, n1, n2);
		end

		function plotElem(al, varargin)
			plotElem@Elem(al, 'Color', [0 0.4470 0.7410], 'LineWidth', 1, varargin{:});
		end
	end
end