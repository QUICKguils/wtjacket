function modeling(sdiv, opts)
% MODELING  Model of the wt jacket, using 3D beam elements.
%
% Arguments:
%	sdiv (int)
%	  Number of subdivisions in the bare structure.
%	opts (char {'p', 'w'})
%	  'p' -> Enable plots creation.
%	  'w' -> Write plotting data in external file.
% Save:
%	SS (struct) with fields:
%		listNode {1xN Node} -- Cell list of subdivised structure nodes.
%		listElem {1xN Elem} -- Cell list of subdivised structure elements.
%		nbNode   (int)      -- Number    of subdivised structure nodes.
%		nbElem   (int)      -- Number    of subdivised structure elements.
%		nbDOF    (int)      -- Number    of subdivised structure DOFs.

%% TODO

% - Implement the write option.
% - Tidy up shared variables and argument variables.

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and bare structure data.
C  = load(fullfile(file_dir, "../res/constants.mat"));
BS = load(fullfile(file_dir, "../res/bare_struct.mat"));

%% Mains solve

% 1. Subdivised structure

% Lists of subnodes and subelements that are created.
SS = sdiv_struct(sdiv);

% Plotting the subdivised structure, if desired.
if contains(opts, 'p')
	plot_sdiv_struct(SS.listNode, SS.listElem);
end

% 2. K_el and M_el

% List of elementary K and M matrices, in local axes.
[K_el, M_el] = set_el_list(SS.listElem);

% 3. K_es and M_es

[K_es, M_es] = set_es_list(SS.listElem, K_el, M_el);

% 4. K and M

% K_free and M_free.
[K_free, M_free] = set_global_matrices(SS, K_es, M_es);

% Apply boundary conditions to obtain K and M.
% WARN: do not forget lumped nacelle mass.

% 5. Solve the eigenvalue problem.
% HINT: use eigs, with 'smallestabs' option.

% 6. Eigenmodes plot

% 7. Convergence study

% 8. Total mass sanity check

% 9. Save data into sdiv_struct.mat

save(fullfile(file_dir, "../res/sdiv_struct.mat"), "-struct", "SS");

%% 1. Subdivised structure

	function SS = sdiv_struct(sdiv)
		% SDIV_STRUCT  Generate nodes and elements of the subdivised structure.
		%
		% Argument:
		%	sdiv (int) -- Number of subsivisions desired.
		% Return:
		%	SS (struct) with fields:
		%		listNode {1xN Node} -- Cell list of nodes.
		%		listElem {1xN Elem} -- Cell list of elements.
		%		nbNode   (int)      -- Number of nodes.
		%		nbElem   (int)      -- Number of elements.
		%		nbDOF    (int)      -- Number of DOFs.

		% Preallocate the node and element lists of the subdivised structure.
		nbNode   = BS.nbNode + (sdiv-1)*BS.nbElem;
		nbElem   = sdiv*BS.nbElem;
		listNode = cell(1, nbNode);
		listElem = cell(1, nbElem);
		% Prefill the listNode with existing bare nodes.
		listNode(1:BS.nbNode) = BS.listNode;

		for i = 1:BS.nbElem
			% Extract element and extremity nodes (for terseness).
			elem = BS.listElem{i};
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
			% Determine link type of selem.
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
			% Fill the listNode and listElem.
			% TODO: try to find a more readable way of doing this.
			listNode(BS.nbNode + (i-1)*(numel(snode)-2) + (1:(numel(snode)-2))) = snode(2:end-1);
			listElem((i-1)*(numel(selem)) + (1:numel(selem))) = selem;

		end

		% Build return data structure.
		SS.listNode = listNode;
		SS.listElem = listElem;
		SS.nbNode   = nbNode;
		SS.nbElem   = nbElem;
		SS.nbDOF    = nbNode * Node.nbDOF;
	end

	function plot_sdiv_struct(listNode, listElem)
		% PLOT_SDIV_STRUCT  Get an overview of the structure.
		%
		% Arguments:
		%	listNode {1xN Node} -- Cell list of nodes.
		%	listElem {1xN Elem} -- Cell list of elements.

		% Instantiate a figure object.
		figure("WindowStyle", "docked");
		hold on;
		% Plot the nodes.
		for node = listNode(1:end)
			if node{:}.pos(3) > C.f_height(end)  % ignone nacelle
				continue
			end
			node{:}.plotNode()
		end
		% Plot the elements.
		for elem = listElem(1:end)
			if elem{:}.n1.pos(3) > C.f_height(end) ...
					|| elem{:}.n2.pos(3) > C.f_height(end)  % ignore nacelle
				continue
			end
			elem{:}.plotElem()
		end
		% Dress the plot.
		xlabel("X-coord. [m]");
		ylabel("Y-coord. [m]");
		zlabel("Z-coord. [m]");
		title("Subdivised structure");
		axis equal;
		grid;
		view([-0.75, -1, 0.75]);
		hold off;
	end

%% 2. K_el and M_el

	% TODO: double check the K_el and M_el matrices.
	function [K_el, M_el] = set_el_matrices(elem)
		% SET_EL_MATRICES  Set the elementary matrices, local axes.
		%
		% Argument:
		%	elem (Elem) -- Structural element.
		% Returns:
		%	K_el (12x12 double) -- Stiffness elementary matrix, local axes.
		%	M_el (12x12 double) -- Mass      elementary matrix, local axes.
		
		% Pre-terms computation for stiffness matrix.
		k = [
			   elem.E*elem.area / elem.length,      elem.G*elem.Jx  / elem.length,   ... % k(1), k(2)
			 2*elem.E*elem.Iyy  / elem.length,    2*elem.E*elem.Izz / elem.length,   ... % k(3), k(4)
			 4*elem.E*elem.Iyy  / elem.length,    4*elem.E*elem.Izz / elem.length,   ... % k(5), k(6)
			 6*elem.E*elem.Iyy  / elem.length^2,  6*elem.E*elem.Izz / elem.length^2, ... % k(7), k(8)
			12*elem.E*elem.Iyy  / elem.length^3, 12*elem.E*elem.Izz / elem.length^3];    % k(9), k(10)

        % Elementary stiffness matrix in local axes [N/m].
		 K_el = [
			 k(1), 0,      0,     0,     0,     0,     -k(1), 0,      0,     0,     0,         0;
			    0, k(10),  0,     0,     0,     k(8),  0,     -k(10), 0,     0,     0,      k(8);
			    0, 0,      k(9),  0,     -k(7), 0,     0,     0,      -k(9), 0,     -k(7),     0;
			    0, 0,      0,     k(2),  0,     0,     0,     0,      0,     -k(2), 0,         0;
			    0, 0,      -k(7), 0,     k(5),  0,     0,     0,      k(7),  0,     k(3),      0;
			    0, k(8),   0,     0,     0,     k(6),  0,     -k(8),  0,     0,     0,      k(4);
			-k(1), 0,      0,     0,     0,     0,     k(1),  0,      0,     0,     0,         0;
			    0, -k(10), 0,     0,     0,     -k(8), 0,     k(10),  0,     0,     0,     -k(8);
			    0, 0,      -k(9), 0,     k(7),  0,     0,     0,      k(9),  0,     k(7),      0;
			    0, 0,      0,     -k(2), 0,     0,     0,     0,      0,     k(2),  0,         0;
			    0, 0,      -k(7), 0,     k(3),  0,     0,     0,      k(7),  0,     k(5),      0;
			    0, k(8),   0,     0,     0,     k(4),  0,     -k(8),  0,     0,     0,      k(6)];

		% Pre-terms computation for stiffness matrix.
		m = [
			     (elem.d/2)^2 / 3,     (elem.d/2)^2 / 6,   ... % m(1), m(2)
			                1 / 3,                1 / 6,   ... % m(3), m(4)
			              13 / 35,               9 / 70,   ... % m(5), m(6)
			 elem.length*11 / 210, elem.length*13 / 420,   ... % m(7), m(8)
			  elem.length^2 / 105,  elem.length^2 / 140];      % m(9), m(10)

		% Elementary mass matrix in local axes [kg].
		M_el = elem.mass * [
			m(3),     0,     0,    0,      0,      0, m(4),     0,     0,    0,      0,      0;
			   0,  m(5),     0,    0,      0,   m(7),    0,  m(6),     0,    0,      0,  -m(8);
			   0,     0,  m(5),    0,  -m(7),      0,    0,     0,  m(6),    0,   m(8),      0;
			   0,     0,     0, m(1),      0,      0,    0,     0,     0, m(2),      0,      0;
			   0,     0, -m(7),    0,   m(9),      0,    0,     0, -m(8),    0, -m(10),      0;
			   0,  m(7),     0,    0,      0,   m(9),    0,  m(8),     0,    0,      0, -m(10);
			m(4),     0,     0,    0,      0,      0, m(3),     0,     0,    0,      0,      0;
			   0,  m(6),     0,    0,      0,   m(8),    0,  m(5),     0,    0,      0,  -m(7);
			   0,     0,  m(6),    0,  -m(8),      0,    0,     0,  m(5),    0,   m(7),      0;
			   0,     0,     0, m(2),      0,      0,    0,     0,     0, m(1),      0,      0;
			   0,     0,  m(8),    0, -m(10),      0,    0,     0,  m(7),    0,   m(9),      0;
			   0, -m(8),     0,    0,      0, -m(10),    0, -m(7),     0,    0,      0,   m(9)];
	end

	function [K_el, M_el] = set_el_list(listElem)
		% SET_EL_LIST  Set all the elementary matrices, local axes.
		%
		% Argument:
		%	listElem {1xN Elem} -- Cell list of elements.
		% Returns:
		%	K_el {1xN (12x12 double)} -- Cell list of K_el matrices.
		%	M_el {1xN (12x12 double)} -- Cell list of M_el matrices.

		K_el = cell(1, numel(listElem));
		M_el = cell(1, numel(listElem));

		for i = 1:numel(listElem)
			elem = listElem{i};
			[K_el{i}, M_el{i}] = set_el_matrices(elem);
		end
	end

%% 3. K_es and M_es

	function [K_es, M_es] = set_es_matrices(elem, K_el, M_el)
		% SET_ES_MATRICES  Set the elementary matrices, structural axes.
		%
		% Arguments:
		%	elem (Elem)         -- Structural element.
		%	K_el (12x12 double) -- Stiffness elementary matrix, local axes.
		%	M_el (12x12 double) -- Mass      elementary matrix, local axes.
		% Returns:
		%	K_el (12x12 double) -- Stiffness elementary matrix, structural axes.
		%	M_el (12x12 double) -- Mass      elementary matrix, structural axes.

		% Build the local basis.
		%
		% Normalized x-axis directional vector of the local axes.
		ex = elem.dir';
		% WARN: see if left-hand bases are OK.
		% Generate ey and ez by computing the null space of {ex, 0, 0}.
		eyz = null([ex, zeros(3, 1), zeros(3, 1)]);
		ey = eyz(:, 1);
		ez = eyz(:, 2);
		% WARN: verify if I should take the transpose.
		lbasis = [ex, ey, ez];

		% Transformation matrix between local and structural axes.
		% The local basis happens to be the rotation operator, as the chosen
		% structural basis is simply the identity matrix.
		T = kron(eye(4), lbasis);

		% Applying the change of basis to K_el and M_el.
		K_es = T' * K_el * T;
        M_es = T' * M_el * T;
	end

	function [K_es, M_es] = set_es_list(listElem, K_el, M_el)
		% SET_ES_LIST  Set all the K and M elementary matrices, structural axes.
		%
		% Arguments:
		%	listElem {1xN Elem}           -- Cell list of elements.
		%	K_el     {1xN (12x12 double)} -- Cell list of K_el matrices.
		%	M_el     {1xN (12x12 double)} -- Cell list of M_el matrices.
		% Returns:
		%	K_es {1xN (12x12 double)} -- Cell list of K_es matrices.
		%	M_es {1xN (12x12 double)} -- Cell list of M_es matrices.

		K_es = cell(1, numel(listElem));
		M_es = cell(1, numel(listElem));

		for i = 1:numel(listElem)
			[K_es{i}, M_es{i}] = set_es_matrices(listElem{i}, K_el{i}, M_el{i});
		end
	end

%% 4. K and M

	function [K_free, M_free] = set_global_matrices(SS, K_es, M_es)
		% SET_GLOBAL_MATRICES  Set the global matrices of the free structure.
		%
		% Arguments:
		%	SS   (struct)             -- Subdivised structure.
		%	K_es {1xN (12x12 double)} -- Cell list of K_es matrices.
		%	M_es {1xN (12x12 double)} -- Cell list of M_es matrices.
		% Returns:
		%	K_free  (nbDOFxnbDOF double) -- Global free stiffness matrix.
		%	M_free  (nbDOFxnbDOF double) -- Global free mass      matrix.

		K_free = zeros(SS.nbDOF);
		M_free = zeros(SS.nbDOF);

		for i = 1:SS.nbElem
			locel = [SS.listElem{i}.n1.dof, SS.listElem{i}.n2.dof];
			K_free(locel, locel) = K_free(locel, locel) + K_es{i};
			M_free(locel, locel) = M_free(locel, locel) + M_es{i};
		end
	end

	function [K, M] = apply_bc(K_free, M_free)
		% APPLY_BC  Apply the boundary conditions.
		%
		% arguments:
		%	K_free  (nbDOFxnbDOF double) -- Global free stiffness matrix.
		%	M_free  (nbDOFxnbDOF double) -- Global free mass      matrix.
		% Returns:
		%	K (nbDOFxnbDOF double) -- Global stiffness matrix.
		%	M (nbDOFxnbDOF double) -- Global mass      matrix.
	end
end
