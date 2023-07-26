classdef RigidLink < MainLink
	% RIGIDLINK  Represent a rigid link.
	%	Rigid links are derived from main links. The material and
	%	geometric parameters are set in such a way that the stiffness is
	%	noticeably increased, thus behaving like a rigid link.

	methods
		function rl = RigidLink(n1, n2)
			% RIGIDLINK  Construct an instance of RigidLink.

			rl = rl@MainLink(n1, n2);

			rl.rho  = rl.rho  * 1e-4;
			rl.E    = rl.E    * 1e4;
			rl.area = rl.area * 1e-2;
			rl.Iyy  = rl.Iyy  * 1e4;
			rl.Izz  = rl.Izz  * 1e4;
			rl.Jx   = rl.Jx   * 1e4;
		end


		function plotElem(rl)
			plotElem@Elem(rl, Color=[0.4660 0.6740 0.1880], LineWidth=1.5);
		end
	end
end
