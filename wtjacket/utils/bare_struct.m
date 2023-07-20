function bare_struct(opts)
% BARE_STRUCT  Definition of the wind turbine jacket structure.
%
% Argument:
%	opts (char {'p'}) -- Enable plots creation.
% Save:
%	BS (Struct) with fields:
%	  nodeList {1xN Node} -- Cell of bare structure nodes.
%	  elemList {1xN Elem} -- Cell of bare structure elements.

% Reset internal state of classes (counters, etc).
clear Node;
clear Elem;

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Load project constant data.
C  = load(fullfile(file_dir, "../../res/constants.mat"));

%% Nodes

	function frame = elevate(h)
		% ELEVATE  Create a frame of nodes, for the desired altitude.
		%
		% Argument: h       (double)   -- Altitude [m].
		% Return:   elevate {1x4 Node} -- Frame of nodes.
		shift = tand(C.angle) * h;
		frame = {
			Node([          shift,           shift,  h]);
			Node([C.b_width-shift,           shift,  h]);
			Node([C.b_width-shift, C.b_width-shift,  h]);
			Node([          shift, C.b_width-shift,  h])};
	end

% Cell containing the 22 nodes of the bare structure.
nodeList = [
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

%% Bare structure plot

	function plot_bare_struct(nl, el)
		% PLOT_BARE_STRUCT  Get an overview of the bare structure.
		%
		% Arguments:
		%	nl {1xN Node} -- Cell of nodes.
		%	el {1xN Elem} -- Cell of elements.

		% Instantiate a figure object.
		figure("WindowStyle", "docked");
		hold on;
		% Plot the nodes.
		for node = nl(1:end-1)  % ignore nacelle
			node{:}.plotNode()
		end
		% Plot the elements.
		for elem = el(1:end-1)  % ignore nacelle
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
	plot_bare_struct(nodeList, elemList);
end

%% Save data into bare_struct.mat

% Gather relevant data to save.
BS.elemList = elemList;
BS.nodeList = nodeList;

% Save data in bare_struct.mat, which lies in the /res directory.
save(fullfile(file_dir, "../../res/bare_struct.mat"), "-struct", "BS");

end
