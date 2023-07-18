classdef Node
	% NODE  Represent a 3D structural node.

	properties
		pos   (1, 3) double {mustBeReal}
		label (1, 1) double {mustBeInteger}
	end

	methods (Static)
		function label = setLabel()
			% SETLABEL  Set the label of the Node instance.
			%   This closure keeps track of the number of Node instances
			%   through a persistent variable.
			persistent cnt;
			if isempty(cnt)
				cnt = 0;
			end
			cnt = cnt + 1;
			label = cnt;
		end
	end

	methods
		function node = Node(pos)
			% NODE  Construct an instance of Node.
			node.pos = pos;
			node.label = node.setLabel;
		end

		function plotNode(node)
			plot3(node.pos(1), node.pos(2), node.pos(3), ...
				Color=[0.4940 0.1840 0.5560], Marker="*", MarkerSize=8, LineWidth=1.5);
		end
	end
end