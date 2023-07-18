function BareStruct(opts)
% BASESTRUCT  Definition of the wind turbine jacket structure.
%
% Argument:
%	opts: char {'p'}, optional. Default is 'p'.
%	  'p' -> Enable plots creation.

% Reset internal state of classes (counters, etc).
clear Node;
clear Elem;

% Option defaults: generate the plots.
if ~nargin
	opts = 'p';
end

%% Nodes

% X and Y coord of the nodes at the base of the structure.
BasePos = [
	0, 0;
	5, 0;
	5, 5;
	0, 5];

	function Frame = elevate(bp, h)
		% ELEVATE  Create the 4 nodes of a frame at altitude h [m].

		% Imposed structure angle [°].
		angle = 3;

		shift = tand(angle) * h;
		Frame = [
			Node([bp(1, 1)+shift, bp(1, 2)+shift, h]);
			Node([bp(2, 1)-shift, bp(2, 2)+shift, h]);
			Node([bp(3, 1)-shift, bp(3, 2)-shift, h]);
			Node([bp(4, 1)+shift, bp(4, 2)-shift, h])];
	end

% Array containing the 20 nodes of the bare structure.
NodeList = [
	% The four horizontal frames.
	elevate(BasePos,  0);
	elevate(BasePos,  1);
	elevate(BasePos,  9);
	elevate(BasePos, 17);
	elevate(BasePos, 25);
	% The two nodes to attach the nacelle.
	Node([5/2, 5/2, 25]);
	Node([5/2, 5/2, 80])];

% Organize the node coordinates in a 2D array.
nodeList = reshape([NodeList.pos], 3, [])';

%% Elements

ElemList = {
	% First main leg.
	MainLink( NodeList(1),  NodeList(5) );
	MainLink( NodeList(5),  NodeList(9) );
	MainLink( NodeList(9),  NodeList(13));
	MainLink( NodeList(13), NodeList(17));
	% Second main leg.
	MainLink( NodeList(2),  NodeList(6) );
	MainLink( NodeList(6),  NodeList(10));
	MainLink( NodeList(10), NodeList(14));
	MainLink( NodeList(14), NodeList(18));
	% Third main leg.
	MainLink( NodeList(3),  NodeList(7) );
	MainLink( NodeList(7),  NodeList(11));
	MainLink( NodeList(11), NodeList(15));
	MainLink( NodeList(15), NodeList(19));
	% Fourth main leg.
	MainLink( NodeList(4),  NodeList(8) );
	MainLink( NodeList(8),  NodeList(12));
	MainLink( NodeList(12), NodeList(16));
	MainLink( NodeList(16), NodeList(20));
	% First frame.
	AuxLink(  NodeList(5),  NodeList(6));
	AuxLink(  NodeList(6),  NodeList(7));
	AuxLink(  NodeList(7),  NodeList(8));
	AuxLink(  NodeList(8),  NodeList(5));
	% Second frame.
	AuxLink(  NodeList(9),  NodeList(10));
	AuxLink(  NodeList(10), NodeList(11));
	AuxLink(  NodeList(11), NodeList(12));
	AuxLink(  NodeList(12), NodeList(9) );
	% Third frame.
	AuxLink(  NodeList(13), NodeList(14));
	AuxLink(  NodeList(14), NodeList(15));
	AuxLink(  NodeList(15), NodeList(16));
	AuxLink(  NodeList(16), NodeList(13));
	% Third frame.
	AuxLink(  NodeList(17), NodeList(18));
	AuxLink(  NodeList(18), NodeList(19));
	AuxLink(  NodeList(19), NodeList(20));
	AuxLink(  NodeList(20), NodeList(17));
	% First stage.
	AuxLink(  NodeList(6),  NodeList(9) );
	AuxLink(  NodeList(6),  NodeList(11));
	AuxLink(  NodeList(8),  NodeList(9) );
	AuxLink(  NodeList(8),  NodeList(11));
	% Second stage.
	AuxLink(  NodeList(9),  NodeList(14));
	AuxLink(  NodeList(9),  NodeList(16));
	AuxLink(  NodeList(11), NodeList(14));
	AuxLink(  NodeList(11), NodeList(16));
	% Third stage.
	AuxLink(  NodeList(14), NodeList(17));
	AuxLink(  NodeList(14), NodeList(19));
	AuxLink(  NodeList(16), NodeList(17));
	AuxLink(  NodeList(16), NodeList(19));
	% Nacelle support.
	RigidLink(NodeList(17), NodeList(21));
	RigidLink(NodeList(18), NodeList(21));
	RigidLink(NodeList(19), NodeList(21));
	RigidLink(NodeList(20), NodeList(21));
	RigidLink(NodeList(21), NodeList(22));
	};


	function plotBareStructure(el, nl)
		% PLOTELEMLIST  Get an overview of the bare structure.
		figure("WindowStyle", "docked");
		hold on;
		for i = 1:numel(el)-1  % neglect nacelle beam
			el{i}.plotElem
		end
		for i = 1:numel(nl)-1  % neglect nacelle
			nl(i).plotNode
		end
		axis equal; grid; view([-0.75, -1, 0.75]);
		hold off;
	end

if contains(opts, 'p')
	plotBareStructure(ElemList, NodeList);
end

%% Save data into BareStruct.mat

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Save data in constants.mat, which lies in the root directory.
save(fullfile(file_dir, "../../res/BareStruct.mat"), "NodeList", "nodeList", "ElemList");

end