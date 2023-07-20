classdef Node
	% NODE  Represent a 3D structural node.

	properties (Constant)
		ndof = 6;  % Number of DOFs of a node.
	end

	properties
		pos   (1, 3) double {mustBeReal}
		label (1, 1) double {mustBeInteger}
		dof   (1, 6) double {mustBeInteger}
	end

	methods (Static)
		function label = setLabel()
			% SETLABEL  Set the label of the Node instance.
			%	This method keeps track of the label number assignment
			%	through a persistent variable.
			persistent label_cnt;
			if isempty(label_cnt)
				label_cnt = 0;
			end
			label_cnt = label_cnt + 1;
			label = label_cnt;
		end

		function dof = setDof()
			% SETDOF  Set the DOFs of the Node instance.
			%	This method keeps track of the DOFs number assigment
			%	though a persistent variable.
			persistent dof_cnt;
			if isempty(dof_cnt)
				dof_cnt = 0;
			end
			dof = dof_cnt + (1:Node.ndof);
			dof_cnt = dof_cnt + Node.ndof;
		end
	end

	methods
		function node = Node(pos)
			node.pos   = pos;
			node.label = node.setLabel;
			node.dof   = node.setDof;
		end

		function plotNode(node)
			plot3(node.pos(1), node.pos(2), node.pos(3), ...
				Color=[0.4940 0.1840 0.5560], Marker="*", MarkerSize=8, LineWidth=1.5);
		end
	end
end