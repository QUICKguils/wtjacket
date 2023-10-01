function modeling(sdiv, plt)
% MODELING  Model of the wt jacket, using 3D beam elements.
%
% Arguments:
%	sdiv (int)            -- Number of subdivisions in the bare structure.
%	plt  (char {'p', ''}) -- 'p' -> Enable plots creation.

prepare_fem_simulation(plt);
run_fem_simulation(sdiv, plt);
end

function prepare_fem_simulation(plt)
% PREPARE_FEM_SIMULATION  Set the FEM simulation initial state.
%
%	Argument:
%	  plt (char {'p', ''}) -- 'p' -> Enable plots creation.

% Reset class internal states, close previous plots.
clear Node Elem
close all;

% Initialize MAT file.
constants();
bare_struct(plt);
end

function run_fem_simulation(sdiv, plt)
% RUN_FEM_SIMULATION  Run a finite element method simulation.
%
% Arguments:
%	sdiv (int)            -- Number of subdivisions in the bare structure.
%	plt  (char {'p', ''}) -- 'p' -> Enable plots creation.
% Save:
%	SS (struct) -- Subdivised structure, with fields:
%		listNode {1xN Node} -- Cell list of nodes.
%		listElem {1xN Elem} -- Cell list of elements.
%		nbNode   (int)      -- Number    of nodes.
%		nbElem   (int)      -- Number    of elements.
%		nbDOF    (int)      -- Number    of DOFs.
%	SOL (struct) -- Solution of the vibration problem, with fields:
%		modes  (nbDOF x nbMode double) -- Modal displacement vectors.
%		freqs  (1 x nbMode)            -- Natural frequencies, in Hertz.
%		nbMode (int)                   -- Number of computed modes.

%TODO:
% - Elementary matrices should be appended to the underlying element,
%   and not defined in detached cell arrays that rely on synced indexing.
% - Tidy up shared variables and argument variables. (for ex; why BS is
%   shared and not SS ?)

%% Imports

% Directory where the present file lies.
file_dir = fileparts(mfilename("fullpath"));

% Constants and bare structure data.
C  = load(fullfile(file_dir, "../res/constants.mat"));
BS = load(fullfile(file_dir, "../res/bare_struct.mat"));

%% Mains solve

% 1. Subdivised structure

SS = sdiv_struct(sdiv);

% Plotting the subdivised structure, if desired.
if contains(plt, 'p')
	plot_sdiv_struct(SS.listNode, SS.listElem);
end

% 2. K_el and M_el

[K_el, M_el] = set_el_list(SS.listElem);

% 3. K_es and M_es

[K_es, M_es] = set_es_list(SS.listElem, K_el, M_el);

% 4. K and M

% K_free and M_free.
[K_free, M_free] = set_global_matrices(SS, K_es, M_es);

% Gather the constrained DOFs in a bool indexing array.
maskCstr = create_cstr_mask(SS);

% Apply boundary conditions.
[K, M] = apply_bc(K_free, M_free, maskCstr);

% 5. Eigenvalue problem solving

SOL = get_solution(K, M, SS, maskCstr);

% 6. Eigenmodes plot

if contains(plt, 'p')
	plot_vibration_mode(SS, SOL);
end

% 7. Total mass and sanity checks

SOL.mass_rbm = rbm_checks(K_free, M_free, SS.nbNode);

% 8. Data saving

% TODO: maybe remove el and es matrices, when fully debugged.
% They are not useful for the following, and take a rather alrge amount
% of space.
KM.K_el = K_el;
KM.M_el = M_el;
KM.K_es = K_es;
KM.M_es = M_es;
KM.K_free = K_free;
KM.M_free = M_free;
KM.K = K;
KM.M = M;

save(fullfile(file_dir, "../res/sdiv_struct.mat"), "-struct", "SS");
save(fullfile(file_dir, "../res/modeling_sol.mat"), "-struct", "SOL");
save(fullfile(file_dir, "../res/modeling_mat.mat"), "-struct", "KM");

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
			if node{:}.pos(3) > C.frame_height(end)  % ignone nacelle
				continue
			end
			node{:}.plotNode()
		end
		% Plot the elements.
		for elem = listElem(1:end)
			if elem{:}.n1.pos(3) > C.frame_height(end) ...
					|| elem{:}.n2.pos(3) > C.frame_height(end)  % ignore nacelle
				continue
			end
			elem{:}.plotElem()
		end
		% Dress the plot.
		xlabel("X/m");
		ylabel("Y/m");
		zlabel("Z/m");
		title("Subdivised structure");
		axis equal;
		grid;
		view(-35, 40);
		hold off;
	end

%% 2. K_el and M_el

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

        % Elementary stiffness matrix in local axes.
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

		% Elementary mass matrix in local axes.
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
			[K_el{i}, M_el{i}] = set_el_matrices(listElem{i});
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
		%	K_es (12x12 double) -- Stiffness elementary matrix, structural axes.
		%	M_es (12x12 double) -- Mass      elementary matrix, structural axes.

		% Build the local basis.
		%
		% Normalized x-axis directional vector of the local axes.
		ex = elem.dir';
		% Generate ey and ez by computing the null space of {ex, 0, 0}.
		nullspace = null(ex');
		ey = nullspace(:, 1);
		ez = cross(ex, ey);  % Not nullspace(:,2), to ensure right-handedness.
		lbasis = [ex, ey, ez];
		
		% Transformation matrix between local and structural axes.
		% The transpose local basis happens to be the rotation operator,
		% as the chosen structural basis is simply the identity matrix.
		T = kron(eye(4), lbasis');

		% Applying the change of basis to K_el and M_el.
		K_es = T' * K_el * T;
        M_es = T' * M_el * T;

		check_sym(K_es);
		check_sym(M_es);
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

		% Add the nacelle inertia.
		% TODO: not super robust to hardcode the nacelle index.
		tr_dofs  = BS.listNode{end}.dof(1:3);
		rot_dofs = BS.listNode{end}.dof(4:6);
		tr_idx  = sub2ind(size(M_free), tr_dofs, tr_dofs);
		rot_idx = sub2ind(size(M_free), rot_dofs, rot_dofs);
		M_free(tr_idx)  = M_free(tr_idx) + C.nacelle_mass;
		M_free(rot_idx) = M_free(rot_idx) + C.nacelle_inertia;

		check_sym(K_free);
		check_sym(M_free);
	end

	function maskCstr = create_cstr_mask(SS)
		% CREATE_CSTR_MASK  Create a bool array that index on constrained DOFs.
		%
		% Argument:
		%	SS (struct) -- Subdivised structure.
		% Return:
		%	maskCstr (1 x nbDOF bool) -- Index on constrained DOFs.

		maskCstr = false(1, SS.nbDOF);

		for node = SS.listNode
			if node{:}.cstr == "clamped"
				maskCstr(node{:}.dof) = true;
			end
		end
	end

	function [K, M] = apply_bc(K_free, M_free, maskCstr)
		% APPLY_BC  Apply the boundary conditions.
		%
		% arguments:
		%	K_free   (nbDOF x nbDOF double) -- Global free stiffness matrix.
		%	M_free   (nbDOF x nbDOF double) -- Global free mass matrix.
		%	maskCstr (1 x nbDOF bool)       -- Index on constrained DOFs.
		% Returns:
		%	K (nbDOF x nbDOF double) -- Global stiffness matrix.
		%	M (nbDOF x nbDOF double) -- Global mass matrix.

		% Supress rows and columns corresponding to clamped node DOFs.
		K_free(:, maskCstr) = [];
		M_free(:, maskCstr) = [];
		K_free(maskCstr, :) = [];
		M_free(maskCstr, :) = [];

		K = K_free;
		M = M_free;
	end

%% 5. Solve the eigenvalue problem

	function SOL = get_solution(K, M, SS, maskCstr, nbMode)
		% GET_SOLUTION  Get the natural frequencies and corresponding modes.
		%
		% Arguments:
		%	K        (nbDOF x nbDOF double) -- Global stiffness matrix.
		%	M        (nbDOF x nbDOF double) -- Global mass matrix.
		%	SS       (struct)               -- Subdivised structure.
		%	maskCstr (1 x nbDOF bool)       -- Index on constrained DOFs.
		%	nbMode   (int, default: 8)      -- Number of first modes desired.
		% Returns:
		%	SOL (struct) -- Solution of the vibration problem, with fields:
		%	  modes  (nbDOF x nbMode double) -- Modal displacement vectors.
		%	  freqs  (1 x nbMode)            -- Natural frequencies, in Hertz.
		%	  nbMode (int)                   -- Number of computed modes.

		% By default, the eight first modes are desired.
		if nargin == 4
			nbMode = 8;
		end

		% Find the first eigvecs and eigvals.
		sigma = 'smallestabs';  % Could be advised to use a scalar value.
		[eigvecs, eigvals] = eigs(K, M, nbMode, sigma);

		% Extract the natural frequencies.
		freqs = sqrt(diag(eigvals)) / (2*pi);

		% Build modal displacements: add clamped nodes.
		modes = zeros(SS.nbDOF, nbMode);
		modes(~maskCstr, :) = eigvecs;

		SOL.freqs  = freqs;
		SOL.modes  = modes;
		SOL.nbMode = nbMode;
	end

%% 6. Eigenmodes plot

	% TODO: Use proper shape function to interpolate the displacement field.
	function plot_vibration_mode(SS, SOL)
		% PLOT_VIBRATION_MODE  Get an overview of the vibration modes.
		%
		% Arguments:
		%	SS  (struct) -- Subdivised structure.
		%	SOL (struct) -- Solution of the vibration problem.

		figure("WindowStyle", "docked");

		% Scale factor.
		ref_length  = C.frame_height(end);
		filter_tr   = reshape((1:3)' + Node.nbDOF * ((1:SS.nbNode)-1), 1, []);
		max_def     = @(idxMode) max(abs(SOL.modes(filter_tr, idxMode)));
		percent_def = 15;  % This gives a readable deformation.
		get_scale   = @(idxMode) ref_length/max_def(idxMode) * percent_def*1e-2;

		for i = 1:SOL.nbMode
			subplot(2, 4, i);
			hold on;

			scale = get_scale(i);
			for elem = SS.listElem
				elem{:}.plotElem();  % pass `'Color', [0, 0, 0, 0.2]` for better clarity.
				x = [elem{:}.n1.pos(1) + scale * SOL.modes(elem{:}.n1.dof(1), i),...
					 elem{:}.n2.pos(1) + scale * SOL.modes(elem{:}.n2.dof(1), i)];
				y = [elem{:}.n1.pos(2) + scale * SOL.modes(elem{:}.n1.dof(2), i), ...
					 elem{:}.n2.pos(2) + scale * SOL.modes(elem{:}.n2.dof(2), i)];
				z = [elem{:}.n1.pos(3) + scale * SOL.modes(elem{:}.n1.dof(3), i), ...
					 elem{:}.n2.pos(3) + scale * SOL.modes(elem{:}.n2.dof(3), i)];
				plot3(x, y, z, Color=[0.9290 0.6940 0.1250], LineWidth=2);
			end

			xlabel("X/m");
			ylabel("Y/m");
			zlabel("Z/m");
			title(['f = ', num2str(SOL.freqs(i)), ' Hz']);
			axis equal;
			grid;
			view(-35, 50);
			hold off;
		end
	end

%% 7. Total mass and sanity checks

	function varargout = rbm_checks(M_free, K_free, nbNode)
		% RBM_CHECKS  Perform sanity checks based on rbm movement.
		%
		% This function uses a translation along the X-axis as a rigid
		% body mode, to check that:
		%   - the total mass of the structure is the same as the one calculated
		%     with this RBM,
		%   - the generalized forces required to maintain the structure in an
		%     overall translational configuration is null.
		%
		% Arguments:
		%	K_free (nbDOF x nbDOF double) -- Global free stiffness matrix.
		%	M_free (nbDOF x nbDOF double) -- Global free mass      matrix.
		%	nbNode (int)                  -- Number of structural nodes.
		% Return:
		%	mass_rbm (double, optional) -- Mass calculated from RBM [kg].

		% Translation of 1m along the X-axis.
		u_rbm = repmat([1, 0, 0, 0, 0, 0]', nbNode, 1);

		% Mass calculated from this translation.
		mass_rbm = u_rbm' * M_free * u_rbm;
		% Forces calculated from this translation.
		g_rbm = K_free * u_rbm;

		% Sanity checks.
		% FIX: this is wrong
		if abs((mass_rbm-BS.mass)/BS.mass) > 1e-2
			warning('wtjacket:WrongRbmMass', ...
				"RBM in translation yields to a wrong total mass calculation.");
		end
		if ~all(g_rbm < 1e-2)
			warning('wtjacket:WrongRbmForces', ...
				"Generalized forces required to translate along a RBM should be null.");
		end
		disp(mass_rbm);

		% Return mass_rbm, if wanted.
		varargout = cell(nargout, 1);
		for k = 1:nargout
			varargout{k} = mass_rbm;
		end
	end
end
