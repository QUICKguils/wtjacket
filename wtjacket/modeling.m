function [BS, SS, KM, SOL] = modeling(C, sdiv, nMode, opts)
% MODELING  Model of the wt jacket, using 3D beam elements.
%
% Arguments:
%	C     (struct)   -- Constant project quantities.
%	sdiv  (int)      -- Number of subdivisions in the bare structure.
%	nMode (int)      -- Number of first mode computed.
%	opts  (1xN char) -- Options.
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
% Return:
%	BS (struct) -- Bare structure, with fields:
%		nodeList {1xN Node}             -- Cell list of nodes.
%		elemList {1xN Elem}             -- Cell list of elements.
%		cmList   {1xN ConcentratedMass} -- Cell list of concentrated masses.
%		nNode    (int)                  -- Number of nodes.
%		nElem    (int)                  -- Number of elements.
%		nDof     (int)                  -- Number of DOFs.
%		mass     (int)                  -- Total mass of the structure [kg].
%	SS (struct) -- Subdivised structure, with fields:
%	  nodeList {1xN Node} -- Cell list of nodes.
%	  elemList {1xN Elem} -- Cell list of elements.
%	  nNode   (int)       -- Number of nodes.
%	  nElem   (int)       -- Number of elements.
%	  nDof    (int)       -- Number of DOFs.
%	KM (struct) -- Global structural matrices, with fields:
%	  KM.K_free (nDof x nDof) -- Global siffness matrix, without constraints.
%	  KM.M_free (nDof x nDof) -- Global mass     matrix, without constraints.
%	  KM.K      (N x N)       -- Global siffness matrix, with    constraints.
%	  KM.M      (N x N)       -- Global mass     matrix, with    constraints.
%	SOL (struct) -- Solution of the vibration problem, with fields:
%	  modes       (nDof x nbMode double) -- Modal displacement vectors.
%	  freqs       (1 x nbMode)           -- Natural frequencies, in Hertz.
%	  nMode       (int)                  -- Number of computed modes.
%	  massFromRbm (double)               -- Mass calculated from RBM [kg].

% Reset class internal states, close previous plots.
clear Node Elem;
close all;

% 1. Subdivised structure

BS = bareStructure(C, opts);
SS = subdivisedStructure(BS, sdiv);

if contains(opts, 'p')
	plotSubdivisedStructure(SS.nodeList, SS.elemList, C.FRAME_HEIGHT(end));
end

% 2. K and M

[KM.K_free, KM.M_free] = buildGlobalMatrices(SS.elemList, BS.cmList, SS.nDof);

cstrMask = buildCstrMask(SS.nodeList, SS.nDof);

[KM.K, KM.M] = applyBoundaryConditions(KM.K_free, KM.M_free, cstrMask);

% 3. Eigenvalue problem solving

[SOL.frequencies, SOL.modes] = solveEigenvalueProblem(KM.K, KM.M, SS.nDof, cstrMask, nMode, 's');
SOL.nMode = nMode;  % Just to keep track of this user's input.

% 4. Eigenmodes plot

if contains(opts, 'p')
	plotVibrationMode(SS.elemList, SS.nNode, SOL, C.FRAME_HEIGHT(end));
end

% 5. Total mass and sanity checks

SOL.massFromRbm = checkRbm(KM.K_free, KM.M_free, SS.nNode, BS.mass);

end

%% 1. Subdivised structure

function SS = subdivisedStructure(BS, sdiv)
% SUBDIVISEDSTRUCTURE  Generate nodes and elements of the subdivised structure.
%
% Argument:
%	sdiv (int)    -- Number of subsivisions desired.
%	BS   (struct) -- Bare structure.
% Return:
%	SS (struct) -- Subdivised structure, with fields:
%	  nodeList {1xN Node} -- Cell list of nodes.
%	  elemList {1xN Elem} -- Cell list of elements.
%	  nNode   (int)       -- Number of nodes.
%	  nElem   (int)       -- Number of elements.
%	  nDof    (int)       -- Number of DOFs.

% Preallocate the node and element lists of the subdivised structure.
nNode   = BS.nNode + (sdiv-1) * BS.nElem;
nElem   = sdiv * BS.nElem;
nodeList = cell(1, nNode);
elemList = cell(1, nElem);
% Prefill the nodeList with existing bare nodes.
nodeList(1:BS.nNode) = BS.nodeList;

for i = 1:BS.nElem
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
	% Determine link type of selem.
	switch class(elem)
		case 'MainLink'
			subLink = @MainLink;
		case 'AuxLink'
			subLink = @AuxLink;
		case 'RigidLink'
			subLink = @RigidLink;
		otherwise
			error("Cannot determine the element type.");
	end
	% Create subelements.
	for j = 1:numel(selem)
		selem{j} = subLink(snode{j}, snode{j+1});
	end
	% Fill the nodeList and elemList.
	% TODO: try to find a more readable way of doing this.
	nodeList(BS.nNode + (i-1)*(numel(snode)-2) + (1:(numel(snode)-2))) = snode(2:end-1);
	elemList((i-1)*(numel(selem)) + (1:numel(selem))) = selem;

end

% Build return data structure.
SS.nodeList = nodeList;
SS.elemList = elemList;
SS.nNode    = nNode;
SS.nElem    = nElem;
SS.nDof     = nNode * Node.nDof;
end

function plotSubdivisedStructure(nodeList, elemList, maxHeight)
% PLOTSUBDIVISEDSTRUCTURE  Get an overview of the structure.
%
% Arguments:
%	nodeList  {1xN Node} -- Cell list of nodes.
%	elemList  {1xN Elem} -- Cell list of elements.
%	maxHeight (double)   -- Maximum plotting height [m].

% Instantiate a figure object.
figure("WindowStyle", "docked");
hold on;
% Plot the nodes.
for node = nodeList(1:end)
	if node{:}.pos(3) > maxHeight  % ignore mast and nacelle
		continue
	end
	node{:}.plotNode()
end
% Plot the elements.
for elem = elemList(1:end)
	if elem{:}.n1.pos(3) > maxHeight ...
			|| elem{:}.n2.pos(3) > maxHeight  % ignore mast and nacelle
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

%% 2. K and M

function [K_free, M_free] = buildGlobalMatrices(elemList, cmList, nDof)
% BUILDGLOBALMATRICES  Set the global matrices of the free structure.
%
% Arguments:
%	elemList {1xN Elem}             -- Cell list of elements.
%	cmList   {1xN ConcentratedMass} -- Cell list of concentrated masses.
%	nDof     (int)                  -- Number of structural DOFs.
% Returns:
%	K_free  (nDofxnDof double) -- Global free stiffness matrix.
%	M_free  (nDofxnDof double) -- Global free mass      matrix.

K_free = zeros(nDof);
M_free = zeros(nDof);

% Add the beams.
for elem = elemList
	locel = [elem{:}.n1.dof, elem{:}.n2.dof];
	K_free(locel, locel) = K_free(locel, locel) + elem{:}.K_es;
	M_free(locel, locel) = M_free(locel, locel) + elem{:}.M_es;
end

% Add the concentrated masses.
for cm = cmList
	ind  = sub2ind(size(M_free), cm{:}.node.dof, cm{:}.node.dof);
	M_free(ind)  = M_free(ind) + [repmat(cm{:}.mass, 1, 3), cm{:}.Jx, cm{:}.Iyy, cm{:}.Izz];
end

allclose(K_free, K_free');
allclose(M_free, M_free');
end

function cstrMask = buildCstrMask(nodeList, nDof)
% BUILDCSTRMASK  Create a bool array that index on constrained DOFs.
%
% Argument:
%	nodeList {1xN Node} -- Cell list of nodes.
%	nDof     (int)      -- Number of structural DOFs.
% Return:
%	cstrMask (1 x nDof bool) -- Index on constrained DOFs.

cstrMask = false(1, nDof);

for node = nodeList
	if node{:}.cstr == "clamped"
		cstrMask(node{:}.dof) = true;
	end
end
end

function [K, M] = applyBoundaryConditions(K_free, M_free, cstrMask)
% APPLYBOUNDARYCONDITIONS  Apply the boundary conditions.
%
% arguments:
%	K_free   (nDof x nDof double) -- Global free stiffness matrix.
%	M_free   (nDof x nDof double) -- Global free mass matrix.
%	cstrMask (1 x nDof bool)      -- Index on constrained DOFs.
% Returns:
%	K (nDof x nDof double) -- Global stiffness matrix.
%	M (nDof x nDof double) -- Global mass matrix.

K = K_free;
M = M_free;

% Supress rows and columns corresponding to clamped node DOFs.
K(:, cstrMask) = [];
M(:, cstrMask) = [];
K(cstrMask, :) = [];
M(cstrMask, :) = [];
end

%% 3. Solve the eigenvalue problem

function [frequencies, modes] = solveEigenvalueProblem(K, M, nDof, cstrMask, nMode, solver)
% SOLVEEIGENVALUEPROBLEM  Get the natural frequencies and corresponding modes.
%
% Arguments:
%	K        (nDof x nDof double) -- Global stiffness matrix.
%	M        (nDof x nDof double) -- Global mass matrix.
%	nDof     (int)                -- Number of structural DOFs.
%	cstrMask (1 x nDof bool)      -- Index on constrained DOFs.
%	nMode    (int)                -- Number of first modes desired.
%	solver   (char)               -- Type of solver.
%	  'f' -> Use the eig  solver (better for small full matrices).
%	  's' -> Use the eigs solver (better for large sparse matrices).
% Returns:
%	frequencies (1 x nMode double)    -- Natural frequencies, in Hertz.
%	modes       (nDof x nMode double) -- Modal displacement vectors.

switch solver
	case 's'
		[frequencies, modes] = solveEigs();
	case 'f'
		[frequencies, modes] = solveEig();
	otherwise
		error("Cannot determine which solver to use. Specify either 'f' or 's'.");
end

	function [frequencies, modes] = solveEigs()
		% Solve the eigenvalue problem.
		sigma = 'smallestabs';  % Could be advised to use a scalar value.
		% TODO: try to make the 'isSymmetricDefinite' option working.
		[eigvecs, eigvals] = eigs(K, M, nMode, sigma);

		% Extract the natural frequencies.
		frequencies = sqrt(diag(eigvals)) / (2*pi);

		% Build modal displacements: add clamped nodes.
		modes = zeros(nDof, nMode);
		modes(~cstrMask, :) = eigvecs;
	end

	function [frequencies, modes] = solveEig()
		% Solve the eigenvalue problem.
		[eigvecs, eigvals] = eig(K, M);

		% Sort in ascending order.
		[eigvals, ind] = sort(diag(eigvals));
		eigvecs = eigvecs(:, ind);
		% Extract the first natural frequencies.
		frequencies = sqrt(eigvals) / (2*pi);
		frequencies = frequencies(1:nMode);
		% Extract corresponding eigvecs.
		eigvecs = eigvecs(:, 1:nMode);

		% Build modal displacements: add clamped nodes.
		modes = zeros(nDof, nMode);
		modes(~cstrMask, :) = eigvecs;
	end
end

%% 4. Eigenmodes plot

% TODO: Use proper shape function to interpolate the displacement field.
function plotVibrationMode(elemList, nNode, SOL, referenceLength)
% PLOTVIBRATIONMODE  Get an overview of the vibration modes.
%
% Arguments:
%	elemList        {1xN Node} -- Cell list of elements.
%	nNode           (int)      -- Number of structural nodes.
%	SOL             (struct)   -- Solution of the vibration problem.
%	referenceLength (double)   -- Reference length to scale the mode.

figure("WindowStyle", "docked");

% Scale factor.
filterTranslationDof = reshape((1:3)' + Node.nDof * ((1:nNode)-1), 1, []);
maximumDeformation = @(iMode) max(abs(SOL.modes(filterTranslationDof, iMode)));
percentageDeformation = 15;  % This gives a readable deformation.
computeScale = @(idxMode) referenceLength/maximumDeformation(idxMode) * percentageDeformation*1e-2;

for i = 1:SOL.nMode
	nGraphColumn = 4;
	subplot(ceil(SOL.nMode/nGraphColumn), nGraphColumn, i);
	hold on;

	scale = computeScale(i);
	for elem = elemList
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
	title(['f = ', num2str(SOL.frequencies(i)), ' Hz']);
	axis equal;
	grid;
	view(-35, 50);
	hold off;
end
end

%% 5. Total mass and sanity checks

function varargout = checkRbm(K_free, M_free, nNode, mass)
% CHECKRBM  Perform sanity checks based on rbm movement.
%
% This function uses a translation along the X-axis as a rigid
% body mode, to check that:
%   - the total mass of the structure is the same as the one calculated
%     with this RBM,
%   - the generalized forces required to maintain the structure in an
%     overall translational configuration is null.
%
% Arguments:
%	K_free (nDof x nDof double) -- Global free stiffness matrix.
%	M_free (nDof x nDof double) -- Global free mass      matrix.
%	nNode (int)                 -- Number of structural nodes.
% Return:
%	massFromRbm (double, optional) -- Mass calculated from RBM [kg].

% Rigid translation of 1m along the X-axis.
rbm = repmat([1, 0, 0, 0, 0, 0]', nNode, 1);

% Mass calculated from this translation.
massFromRbm = rbm' * M_free * rbm;
% Forces calculated from this translation.
forceFromRbm = K_free * rbm;

% Sanity checks.
allclose(mass, massFromRbm);
% TODO: implement a proper check for g_rbm

% Return mass_rbm, if wanted.
varargout = cell(nargout, 1);
for k = 1:nargout
	varargout{k} = massFromRbm;
end
end