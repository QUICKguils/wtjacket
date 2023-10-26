classdef Node
	% NODE  Represent a 3D structural node.

	properties (Constant, Hidden)
		nDof     = 6;                                 % Number of DOFs of a node.
		cstrType = categorical(["free", "clamped"]);  % Possible constraints types.
	end

	properties
		pos   (1, 3) double {mustBeReal}     % Node coordinate, in structural axes.
		label (1, 1) double {mustBeInteger}  % Label assigned to the Node.
		dof   (1, 6) double {mustBeInteger}  % DOFs  assigned to the Node.
		cstr  categorical                    % Constraints on the Node.
	end

	methods (Static)
		function label = setLabel()
			% SETLABEL  Set the label of the Node instance.
			%	This method keeps track of the label number assignment
			%	through a persistent variable.
			persistent labelCounter;
			if isempty(labelCounter)
				labelCounter = 0;
			end
			labelCounter = labelCounter + 1;
			label = labelCounter;
		end

		function dof = setDof()
			% SETDOF  Set the DOFs of the Node instance.
			%	This method keeps track of the DOFs number assigment
			%	though a persistent variable.
			persistent dofCounter;
			if isempty(dofCounter)
				dofCounter = 0;
			end
			dof = dofCounter + (1:Node.nDof);
			dofCounter = dofCounter + Node.nDof;
		end
	end

	methods
		function node = Node(pos, cstr)
			node.pos   = pos;
			node.label = node.setLabel;
			node.dof   = node.setDof;
			if nargin == 1
				node.cstr = Node.cstrType(1);  % Default 'free'
			elseif nargin == 2
				node.cstr = Node.cstrType(Node.cstrType == cstr);
				if isempty(node.cstr)
					error("'%s' is not a valid constraint type. Try 'free' or 'clamped'.", cstr);
				end
			end
		end

		function plotNode(node)
			plot3(node.pos(1), node.pos(2), node.pos(3), ...
				Color=[0.4940 0.1840 0.5560], Marker=".", MarkerSize=15);
		end
	end
end
