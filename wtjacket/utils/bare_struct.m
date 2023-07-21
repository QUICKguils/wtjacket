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

	function frame = elevate(h)
		% ELEVATE  Create a frame of nodes, for the desired altitude.
		%
		% Argument: h     (double)   -- Altitude [m].
		% Return:   frame {1x4 Node} -- Frame of nodes.
		shift = tand(C.angle) * h;
		frame = {
			Node([          shift,           shift,  h]);
			Node([C.b_width-shift,           shift,  h]);
			Node([C.b_width-shift, C.b_width-shift,  h]);
			Node([          shift, C.b_width-shift,  h])};
	end

% Cell containing the 22 nodes of the bare structure.
listNode = [
	% The four horizontal frames.
	elevate(C.f_height(1));
	elevate(C.f_height(2));
	elevate(C.f_height(3));
	elevate(C.f_height(4));
	elevate(C.f_height(5));
	% The two nodes to attach the nacelle.
	{Node([C.b_width/2, C.b_width/2, C.f_height(end)])};
	{Node([C.b_width/2, C.b_width/2, C.n_height     ])}]';

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

%% Bare structure plot

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
		xlabel("X-coord. [m]");
		ylabel("Y-coord. [m]");
		zlabel("Z-coord. [m]");
		title("Bare structure");
		axis equal;
		grid;
		view([-0.75, -1, 0.75]);
		hold off;
	end

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

% Save data in bare_struct.mat, which lies in the /res directory.
save(fullfile(file_dir, "../../res/bare_struct.mat"), "-struct", "BS");

end
