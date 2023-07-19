function modeling(sdiv, opts)
% MODELING  Model of the wt jacket, using 3D beam elements.
%
% Arguments:
%	sdiv: int, optional. Default is 4.
%	  Number of subsivisions in the bare structure.
%	opts: char {'p', 'w'}, optional. Default is 'p'.
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Contants and bare structure data.
C  = load(fullfile(file_dir, "../res/constants.mat"));
BS = load(fullfile(file_dir, "../res/bare_struct.mat"));
n_belem = numel(BS.elemList);
n_bnode = numel(BS.nodeList);

%% Options setting

if nargin == 0
	sdiv = 4;
	opts = 'p';
elseif nargin == 1
	opts = 'p';
end

%% Creation of the subdivised structure

% Preallocate the node and element lists of the subdivised structure.
nodeList = cell(1, n_bnode + (sdiv-1)*n_belem);
elemList = cell(1, sdiv*n_belem);
% Prefill the nodeList with existing bare nodes.
nodeList(1: n_bnode) = BS.nodeList;


for i = 1:n_belem
	% Extract element and extremity nodes (for terseness).
	elem = BS.elemList{i};
	n1 = elem.n1;
	n2 = elem.n2;
	% Preallocate new subelements and subnodes.
	selem = cell(1, sdiv);
	snode = cell(1, sdiv+1);
	% Prefill the subnodes.
	snode([1, end]) = {n1, n2};
	% Create subnodes.
	for j = 2:numel(snode)-1
		snode{j} = Node(n1.pos + (n2.pos-n1.pos)./sdiv*(j-1));
	end
	% Determine Link type of selem.
	switch class(elem)
		case 'MainLink'
			subLink = @MainLink;
		case 'AuxLink'
			subLink = @AuxLink;
		case 'RigidLink'
			subLink = @RigidLink;
	end
	% Create subelements.
	for j = 1:numel(selem)
		selem{j} = subLink(snode{j}, snode{j+1});
	end
	% Fill the nodeList and elemList.
	nodeList(n_bnode + (i-1)*(numel(snode)-2) + (1:(numel(snode)-2))) = snode(2:end-1);
	elemList((i-1)*(numel(selem)) + (1:numel(selem))) = selem;

end

	function plot_struct(el, nl)
		% PLOT_STRUCT  Get an overview of the structure.
		figure("WindowStyle", "docked");
		hold on;
		for e = el(1:end)
			if e{:}.n1.pos(3) > C.f_height(end) ...
					|| e{:}.n2.pos(3) > C.f_height(end)  % ignore nacelle
				continue
			end
			e{:}.plotElem()
		end
		for node = nl(1:end)
			if node{:}.pos(3) > C.f_height(end)  % ignone nacelle
				continue
			end
			node{:}.plotNode()
		end
		axis equal; grid; view([-0.75, -1, 0.75]);
		hold off;
	end

if contains(opts, 'p')
	plot_struct(elemList, nodeList);
end

end
