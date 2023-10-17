function BS = bare_structure(C, opts)
% BARE_STRUCTURE  Definition of the wind turbine jacket structure.
%
% Argument:
%	C    (struct)   -- Constant project quantities.
%	opts (1xN char) -- Options.
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
% Return:
%	BS (struct) -- Bare structure, with fields:
%	  nodeList {1xN Node}             -- Cell list of nodes.
%	  elemList {1xN Elem}             -- Cell list of elements.
%	  cmList   {1xN ConcentratedMass} -- Cell list of concentrated masses.
%	  nNode    (int)                  -- Number of nodes.
%	  nElem    (int)                  -- Number of elements.
%	  nDof     (int)                  -- Number of DOFs.
%	  mass     (int)                  -- Total mass of the structure [kg].

%% Nodes

	% TODO: find a more robust way to propagate default argument 'free'.
	function frame = elevate(h, cstrList)
		% ELEVATE  Create a frame of nodes, for the desired altitude.
		%
		% Arguments:
		%	h        (double)        -- Altitude [m].
		%	cstrList (1x4 Node.cstr) -- Array of nodes constraint type (default: 1x4 "free").
		% Return:
		%	frame {1x4 Node} -- Frame of nodes.

		if nargin == 1
			cstrList = repmat("free", 1, 4);
		end

		shift = tand(C.LEG_ANGLE) * h;
		frame = {
			Node([             shift,              shift,  h], cstrList(1));
			Node([C.BASE_WIDTH-shift,              shift,  h], cstrList(2));
			Node([C.BASE_WIDTH-shift, C.BASE_WIDTH-shift,  h], cstrList(3));
			Node([             shift, C.BASE_WIDTH-shift,  h], cstrList(4))};
	end

nodeList = [
	% The four horizontal frames.
	elevate(C.FRAME_HEIGHT(1), repmat("clamped", 1, 4));  % The base is clamped.
	elevate(C.FRAME_HEIGHT(2));
	elevate(C.FRAME_HEIGHT(3));
	elevate(C.FRAME_HEIGHT(4));
	elevate(C.FRAME_HEIGHT(5));
	% The two nodes to attach the nacelle.
	{Node([C.BASE_WIDTH/2, C.BASE_WIDTH/2, C.FRAME_HEIGHT(end)])};
	{Node([C.BASE_WIDTH/2, C.BASE_WIDTH/2, C.NACELLE_HEIGHT   ])}]';

%% Elements

elemList = {
	% First leg.
	MainLink( nodeList{1},  nodeList{5} );
	MainLink( nodeList{5},  nodeList{9} );
	MainLink( nodeList{9},  nodeList{13});
	MainLink( nodeList{13}, nodeList{17});
	% Second leg.
	MainLink( nodeList{2},  nodeList{6} );
	MainLink( nodeList{6},  nodeList{10});
	MainLink( nodeList{10}, nodeList{14});
	MainLink( nodeList{14}, nodeList{18});
	% Third leg.
	MainLink( nodeList{3},  nodeList{7} );
	MainLink( nodeList{7},  nodeList{11});
	MainLink( nodeList{11}, nodeList{15});
	MainLink( nodeList{15}, nodeList{19});
	% Fourth leg.
	MainLink( nodeList{4},  nodeList{8} );
	MainLink( nodeList{8},  nodeList{12});
	MainLink( nodeList{12}, nodeList{16});
	MainLink( nodeList{16}, nodeList{20});
	% First frame.
	AuxLink(  nodeList{5},  nodeList{6});
	AuxLink(  nodeList{6},  nodeList{7});
	AuxLink(  nodeList{7},  nodeList{8});
	AuxLink(  nodeList{8},  nodeList{5});
	% Second frame.
	AuxLink(  nodeList{9},  nodeList{10});
	AuxLink(  nodeList{10}, nodeList{11});
	AuxLink(  nodeList{11}, nodeList{12});
	AuxLink(  nodeList{12}, nodeList{9} );
	% Third frame.
	AuxLink(  nodeList{13}, nodeList{14});
	AuxLink(  nodeList{14}, nodeList{15});
	AuxLink(  nodeList{15}, nodeList{16});
	AuxLink(  nodeList{16}, nodeList{13});
	% Fourth frame.
	AuxLink(  nodeList{17}, nodeList{18});
	AuxLink(  nodeList{18}, nodeList{19});
	AuxLink(  nodeList{19}, nodeList{20});
	AuxLink(  nodeList{20}, nodeList{17});
	% First cell.
	AuxLink(  nodeList{6},  nodeList{9} );
	AuxLink(  nodeList{6},  nodeList{11});
	AuxLink(  nodeList{8},  nodeList{9} );
	AuxLink(  nodeList{8},  nodeList{11});
	% Second cell.
	AuxLink(  nodeList{9},  nodeList{14});
	AuxLink(  nodeList{9},  nodeList{16});
	AuxLink(  nodeList{11}, nodeList{14});
	AuxLink(  nodeList{11}, nodeList{16});
	% Third cell.
	AuxLink(  nodeList{14}, nodeList{17});
	AuxLink(  nodeList{14}, nodeList{19});
	AuxLink(  nodeList{16}, nodeList{17});
	AuxLink(  nodeList{16}, nodeList{19});
	% Nacelle rigid support.
	RigidLink(nodeList{17}, nodeList{21});
	RigidLink(nodeList{18}, nodeList{21});
	RigidLink(nodeList{19}, nodeList{21});
	RigidLink(nodeList{20}, nodeList{21});
	RigidLink(nodeList{21}, nodeList{22});
	}';

%% Concentrated masses

cmList = {ConcentratedMass(nodeList{22})};  % Only one CM: the nacelle.

%% Total mass

% Initialise the total mass of the bare structure [kg].
mass = 0;
% Account for the beams.
for elem = elemList
	mass = mass + elem{:}.mass;
end
% Account for the concentrated masses.
for cm = cmList
	mass = mass + cm{:}.mass;
end

%% Bare structure plot

if contains(opts, 'p')
	plot_bare_structure(nodeList, elemList);
end

%% Build return data structure.

BS.nodeList = nodeList;
BS.elemList = elemList;
BS.cmList   = cmList;
BS.nNode    = numel(nodeList);
BS.nElem    = numel(elemList);
BS.nDof     = BS.nNode * Node.nDof;
BS.mass     = mass;

end

function plot_bare_structure(nodeList, elemList)
% PLOT_BARE_STRUCTURE  Get an overview of the bare structure.
%
% Arguments:
%	nodeList {1xN Node} -- Cell list of nodes.
%	elemList {1xN Elem} -- Cell list of elements.

% Instantiate a figure object.
figure("WindowStyle", "docked");
hold on;
% Plot the nodes.
for node = nodeList(1:end-1)  % ignore nacelle
	node{:}.plotNode()
end
% Plot the elements.
for elem = elemList(1:end-1)  % ignore nacelle
	elem{:}.plotElem()
end
% Dress the plot.
xlabel("X/m");
ylabel("Y/m");
zlabel("Z/m");
title("Bare structure");
axis equal;
grid;
view(-35, 40);
hold off;
end
