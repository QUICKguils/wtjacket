function bare_struct(opts)
% BARE_STRUCT  Definition of the wind turbine jacket structure.
%
% Argument:
%	opts (char {'p', 'w'})
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.
% Save:
%	BS (struct) with fields:
%		listNode {1xN Node} -- Cell list of bare structure nodes.
%		listElem {1xN Elem} -- Cell list of bare structure elements.
%		nbNode   (int)      -- Number    of bare structure nodes.
%		nbElem   (int)      -- Number    of bare structure elements.
%		nbDOF    (int)      -- Number    of bare structure DOFs.

%% TODO

% - Implement the write option.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Load project constant data.
C = load(fullfile(file_dir, "../../res/constants.mat"));

%% Nodes

% TODO: find a more robust way to propagate default argument 'free'.
	function frame = elevate(h, cstrList)
		% ELEVATE  Create a frame of nodes, for the desired altitude.
		%
		% Arguments:
		%	h (double)
		%		Altitude [m].
		%	cstrList (1x4 Node.cstr) -- Optional, default is (1x4 "free").
		%		Array of nodes constraint type.
		% Return:
		%	frame {1x4 Node}
		%		Frame of nodes.

		if nargin == 1
			cstrList = repmat("free", 1, 4);
		end

		shift = tand(C.leg_angle) * h;
		frame = {
			Node([             shift,              shift,  h], cstrList(1));
			Node([C.base_width-shift,              shift,  h], cstrList(2));
			Node([C.base_width-shift, C.base_width-shift,  h], cstrList(3));
			Node([             shift, C.base_width-shift,  h], cstrList(4))};
	end

% Cell containing the 22 nodes of the bare structure.
listNode = [
	% The four horizontal frames.
	elevate(C.frame_height(1), repmat("clamped", 1, 4));  % The base is clamped.
	elevate(C.frame_height(2));
	elevate(C.frame_height(3));
	elevate(C.frame_height(4));
	elevate(C.frame_height(5));
	% The two nodes to attach the nacelle.
	{Node([C.base_width/2, C.base_width/2, C.frame_height(end)])};
	{Node([C.base_width/2, C.base_width/2, C.nacelle_height   ])}]';

%% Elements

listElem = {
	% First leg.
	MainLink( listNode{1},  listNode{5} );
	MainLink( listNode{5},  listNode{9} );
	MainLink( listNode{9},  listNode{13});
	MainLink( listNode{13}, listNode{17});
	% Second leg.
	MainLink( listNode{2},  listNode{6} );
	MainLink( listNode{6},  listNode{10});
	MainLink( listNode{10}, listNode{14});
	MainLink( listNode{14}, listNode{18});
	% Third leg.
	MainLink( listNode{3},  listNode{7} );
	MainLink( listNode{7},  listNode{11});
	MainLink( listNode{11}, listNode{15});
	MainLink( listNode{15}, listNode{19});
	% Fourth leg.
	MainLink( listNode{4},  listNode{8} );
	MainLink( listNode{8},  listNode{12});
	MainLink( listNode{12}, listNode{16});
	MainLink( listNode{16}, listNode{20});
	% First frame.
	AuxLink(  listNode{5},  listNode{6});
	AuxLink(  listNode{6},  listNode{7});
	AuxLink(  listNode{7},  listNode{8});
	AuxLink(  listNode{8},  listNode{5});
	% Second frame.
	AuxLink(  listNode{9},  listNode{10});
	AuxLink(  listNode{10}, listNode{11});
	AuxLink(  listNode{11}, listNode{12});
	AuxLink(  listNode{12}, listNode{9} );
	% Third frame.
	AuxLink(  listNode{13}, listNode{14});
	AuxLink(  listNode{14}, listNode{15});
	AuxLink(  listNode{15}, listNode{16});
	AuxLink(  listNode{16}, listNode{13});
	% Fourth frame.
	AuxLink(  listNode{17}, listNode{18});
	AuxLink(  listNode{18}, listNode{19});
	AuxLink(  listNode{19}, listNode{20});
	AuxLink(  listNode{20}, listNode{17});
	% First cell.
	AuxLink(  listNode{6},  listNode{9} );
	AuxLink(  listNode{6},  listNode{11});
	AuxLink(  listNode{8},  listNode{9} );
	AuxLink(  listNode{8},  listNode{11});
	% Second cell.
	AuxLink(  listNode{9},  listNode{14});
	AuxLink(  listNode{9},  listNode{16});
	AuxLink(  listNode{11}, listNode{14});
	AuxLink(  listNode{11}, listNode{16});
	% Third cell.
	AuxLink(  listNode{14}, listNode{17});
	AuxLink(  listNode{14}, listNode{19});
	AuxLink(  listNode{16}, listNode{17});
	AuxLink(  listNode{16}, listNode{19});
	% Nacelle rigid support.
	RigidLink(listNode{17}, listNode{21});
	RigidLink(listNode{18}, listNode{21});
	RigidLink(listNode{19}, listNode{21});
	RigidLink(listNode{20}, listNode{21});
	RigidLink(listNode{21}, listNode{22});
	}';


%% Total mass

% Initialise the total mass of the bare structure [kg].
mass = 0;
% Account for the beams.
for elem = listElem
	mass = mass + elem{:}.mass;
end
% Account for the nacelle.
mass = mass + C.nacelle_mass;


%% Bare structure plot

if contains(opts, 'p')
	plot_bare_struct(listNode, listElem);
end

%% Save data into bare_struct.mat

% Gather relevant data to save.
BS.listNode = listNode;
BS.listElem = listElem;
BS.nbNode   = numel(listNode);
BS.nbElem   = numel(listElem);
BS.nbDOF    = BS.nbNode * Node.nbDOF;
BS.mass     = mass;

% Save data in bare_struct.mat, which lies in the /res directory.
save(fullfile(file_dir, "../../res/bare_struct.mat"), "-struct", "BS");

end

function plot_bare_struct(listNode, listElem)
% PLOT_BARE_STRUCT  Get an overview of the bare structure.
%
% Arguments:
%	listNode {1xN Node} -- Cell list of nodes.
%	listElem {1xN Elem} -- Cell list of elements.

% Instantiate a figure object.
figure("WindowStyle", "docked");
hold on;
% Plot the nodes.
for node = listNode(1:end-1)  % ignore nacelle
	node{:}.plotNode()
end
% Plot the elements.
for elem = listElem(1:end-1)  % ignore nacelle
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