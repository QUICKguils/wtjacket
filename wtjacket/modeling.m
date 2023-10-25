function [BareStruct, SdivStruct, FemMat, FemSol] = modeling(Cst, sdiv, nMode, opts)
% MODELING  Model of the wt jacket, using 3D beam elements.
%
% Arguments:
%	Cst   (struct)   -- Constant project quantities.
%	sdiv  (int)      -- Number of subdivisions in the bare structure.
%	nMode (int)      -- Number of first mode computed.
%	opts  (1xN char) -- Options.
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
% Return:
%	BareStruct (struct) -- Bare structure, with fields:
%	  nodeList {1 x N Node}             -- Cell list of nodes.
%	  elemList {1 x N Elem}             -- Cell list of elements.
%	  cmList   {1 x N ConcentratedMass} -- Cell list of concentrated masses.
%	  loadList {1 x N Load}             -- Cell list of loads.
%	  nNode    (int)                    -- Number of nodes.
%	  nElem    (int)                    -- Number of elements.
%	  nDof     (int)                    -- Number of DOFs.
%	  mass     (int)                    -- Total mass of the structure [kg].
%	SdivStruct (struct) -- Subdivised structure, with fields:
%	  nodeList {1 x N Node}             -- Cell list of nodes.
%	  elemList {1 x N Elem}             -- Cell list of elements.
%	  loadList {1 x N Load}             -- Cell list of loads.
%	  cmList   {1 x N ConcentratedMass} -- Cell list of concentrated masses.
%	  nNode    (int)                    -- Number of nodes.
%	  nElem    (int)                    -- Number of elements.
%	  nDof     (int)                    -- Number of DOFs.
%	FemMat (struct) -- Global structural matrices, with fields:
%	  K_free   (nDofxnDof double) -- Global siffness matrix, without constraints.
%	  M_free   (nDofxnDof double) -- Global mass     matrix, without constraints.
%	  K        (NxN double)       -- Global siffness matrix, with    constraints.
%	  M        (NxN double)       -- Global mass     matrix, with    constraints.
%	  cstrMask (1xnDof bool)      -- Index on constrained DOFs.
%	FemSol (struct) -- Solution of the FEM simulation, with fields:
%	  frequencyHertz (1 x nMode double)    -- Natural frequencies [Hz].
%	  frequencyRad   (1 x nMode double)    -- Natural frequencies [rad/s].
%	  mode           (nDof x nMode double) -- Modal displacement vectors.
%	  nMode          (int)                 -- Number of computed first modes.
%	  nDof           (int)                 -- Number of DOFs.
%	  massFromRbm    (double)              -- Mass calculated from RBM [kg].

% Reset class internal states.
clear Node Elem;

% 1. Subdivised structure

BareStruct = bare_structure(Cst, opts);
SdivStruct = subdivise_structure(BareStruct, sdiv);

if contains(opts, 'p')
	plot_subdivised_structure(SdivStruct, Cst.FRAME_HEIGHT(end));
end

% 2. K and M

[FemMat.K_free, FemMat.M_free] = build_global_matrices(SdivStruct);

FemMat.cstrMask = build_cstr_mask(SdivStruct);

[FemMat.K, FemMat.M] = apply_boundary_conditions(FemMat);

% 3. Eigenvalue problem solving

FemSol = solve_eigenvalue_problem(SdivStruct, FemMat, nMode, 's');

% 4. Eigenmodes plot

if contains(opts, 'p')
	plot_vibration_mode(SdivStruct, FemSol, Cst.FRAME_HEIGHT(end));
end

% 5. Total mass and sanity checks

FemSol.massFromRbm = check_rbm(FemMat.M_free, SdivStruct.nNode, BareStruct.mass);

end

%% 1. Subdivised structure

function SdivStruct = subdivise_structure(BareStruct, sdiv)
% SUBDIVISE_STRUCTURE  Subdivise the bare structure.
%
% Arguments:
%	BareStruct (struct) -- Bare structure.
%	sdiv       (int)    -- Number of subsivisions desired.
% Return:
%	SdivStruct (struct) -- Subdivised structure, with fields:
%	  nodeList {1 x N Node}             -- Cell list of nodes.
%	  elemList {1 x N Elem}             -- Cell list of elements.
%	  loadList {1 x N Load}             -- Cell list of loads.
%	  cmList   {1 x N ConcentratedMass} -- Cell list of concentrated masses.
%	  nNode    (int)                    -- Number of nodes.
%	  nElem    (int)                    -- Number of elements.
%	  nDof     (int)                    -- Number of DOFs.

% Preallocate the node and element lists of the subdivised structure.
nNode    = BareStruct.nNode + (sdiv-1) * BareStruct.nElem;
nElem    = sdiv * BareStruct.nElem;
nodeList = cell(1, nNode);
elemList = cell(1, nElem);
% Prefill the nodeList with existing bare nodes.
nodeList(1:BareStruct.nNode) = BareStruct.nodeList;

for i = 1:BareStruct.nElem
	% Extract element and extremity nodes (for terseness).
	elem = BareStruct.elemList{i};
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
	nodeList(BareStruct.nNode + (i-1)*(numel(snode)-2) + (1:(numel(snode)-2))) = snode(2:end-1);
	elemList((i-1)*(numel(selem)) + (1:numel(selem))) = selem;

end

% Build return data structure.
SdivStruct.nodeList = nodeList;
SdivStruct.elemList = elemList;
SdivStruct.cmList   = BareStruct.cmList;
SdivStruct.loadList = BareStruct.loadList;
SdivStruct.nNode    = nNode;
SdivStruct.nElem    = nElem;
SdivStruct.nDof     = nNode * Node.nDof;
end

function plot_subdivised_structure(SdivStruct, maxHeight)
% PLOT_SUBDIVISED_STRUCTURE  Get an overview of the structure.
%
% Arguments:
%	SdivStruct (struct) -- Subdivised structure.
%	maxHeight  (double) -- Maximum plotting height [m].

% Instantiate a figure object.
figure("WindowStyle", "docked");
hold on;
% Plot the nodes.
for node = SdivStruct.nodeList(1:end)
	if node{:}.pos(3) > maxHeight  % ignore mast and nacelle
		continue
	end
	node{:}.plotNode()
end
% Plot the elements.
for elem = SdivStruct.elemList(1:end)
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

function [K_free, M_free] = build_global_matrices(SdivStruct)
% BUILD_GLOBAL_MATRICES  Set the global matrices of the free structure.
%
% Argument:
%	SdivStruct (struct) -- Subdivised structure.
% Returns:
%	K_free (nDof x nDof double) -- Global free stiffness matrix.
%	M_free (nDof x nDof double) -- Global free mass      matrix.

K_free = zeros(SdivStruct.nDof);
M_free = zeros(SdivStruct.nDof);

% Add the beams.
for elem = SdivStruct.elemList
	locel = [elem{:}.n1.dof, elem{:}.n2.dof];
	K_free(locel, locel) = K_free(locel, locel) + elem{:}.K_es;
	M_free(locel, locel) = M_free(locel, locel) + elem{:}.M_es;
end

% Add the concentrated masses.
for cm = SdivStruct.cmList
	ind  = sub2ind(size(M_free), cm{:}.node.dof, cm{:}.node.dof);
	M_free(ind)  = M_free(ind) + [repmat(cm{:}.mass, 1, 3), cm{:}.Jx, cm{:}.Iyy, cm{:}.Izz];
end

allclose(K_free, K_free');
allclose(M_free, M_free');
end

function cstrMask = build_cstr_mask(SdivStruct)
% BUILD_CSTR_MASK  Create a bool array that index on constrained DOFs.
%
% Argument:
%	SdivStruct (struct) -- Subdivised structure.
% Return:
%	cstrMask (1 x nDof bool) -- Index on constrained DOFs.

cstrMask = false(1, SdivStruct.nDof);

for node = SdivStruct.nodeList
	if node{:}.cstr == "clamped"
		cstrMask(node{:}.dof) = true;
	end
end
end

function [K, M] = apply_boundary_conditions(FemMat)
% APPLY_BOUNDARY_CONDITIONS  Apply the boundary conditions.
%
% argument:
%	FemMat (struct) -- Global structural matrices.
% Returns:
%	K (N x N double) -- Global stiffness matrix.
%	M (N x N double) -- Global mass matrix.

% Supress rows and columns corresponding to clamped node DOFs.
K = FemMat.K_free(~FemMat.cstrMask, ~FemMat.cstrMask);
M = FemMat.M_free(~FemMat.cstrMask, ~FemMat.cstrMask);
end

%% 3. Solve the eigenvalue problem

function FemSol = solve_eigenvalue_problem(SdivStruct, FemMat, nMode, solver)
% SOLVE_EIGENVALUE_PROBLEM  Get the natural frequencies and corresponding modes.
%
% Arguments:
%	SdivStruct (struct) -- Subdivised structure.
%	FemMat     (struct) -- Global structural matrices.
%	nMode      (int)    -- Number of first modes desired.
%	solver     (char)   -- Type of solver.
%	  'f' -> Use the eig  solver (better for small full matrices).
%	  's' -> Use the eigs solver (better for large sparse matrices).
% Returns:
%	FemSol (struct) -- Solution of the FEM simulation, with fields:
%	  frequencyHertz (1 x nMode double)    -- Natural frequencies, in hertz.
%	  frequencyRad   (1 x nMode double)    -- Natural frequencies, in rad/s.
%	  mode           (nDof x nMode double) -- Modal displacement vectors.
%	  nMode          (int)                 -- Number of computed first modes.
%	  nDof           (int)                 -- Number of DOFs.
switch solver
	case 's'
		[frequencyHertz, frequencyRad, mode] = solve_eigs(nMode);
	case 'f'
		[frequencyHertz, frequencyRad,mode] = solve_eig(nMode);
	otherwise
		error("Cannot determine which solver to use. Specify either 'f' or 's'.");
end

	function [frequencyHetz, frequencyRad, mode] = solve_eigs(nMode)
		% Solve the eigenvalue problem.
		sigma = 'smallestabs';  % Could be advised to use a scalar value.
		% TODO: try to make the 'isSymmetricDefinite' option working.
		[eigvecs, eigvals] = eigs(FemMat.K, FemMat.M, nMode, sigma);

		% Extract the natural frequencies.
		frequencyRad  = sqrt(diag(eigvals));
		frequencyHetz = frequencyRad / (2*pi);

		% Build modal displacements: add clamped nodes.
		mode = zeros(SdivStruct.nDof, nMode);
		mode(~FemMat.cstrMask, :) = eigvecs;
	end

	function [frequencyHertz, frequencyRad, mode] = solve_eig(nMode)
		% Solve the eigenvalue problem.
		[eigvecs, eigvals] = eig(FemMat.K, FemMat.M);

		% Sort in ascending order.
		[eigvals, ind] = sort(diag(eigvals));
		eigvecs = eigvecs(:, ind);
		% Extract the first natural frequencies.
		frequencyRad   = sqrt(eigvals);
		frequencyHertz = frequencyRad / (2*pi);
		frequencyHertz = frequencyHertz(1:nMode);
		% Extract corresponding eigvecs.
		eigvecs = eigvecs(:, 1:nMode);

		% Build modal displacements: add clamped nodes.
		mode = zeros(SdivStruct.nDof, nMode);
		mode(~FemMat.cstrMask, :) = eigvecs;
	end

	FemSol.frequencyHertz = frequencyHertz;
	FemSol.frequencyRad   = frequencyRad;
	FemSol.mode           = mode;
	FemSol.nMode          = nMode;
	FemSol.nDof           = SdivStruct.nDof;
end

%% 4. Eigenmodes plot

% TODO: Use proper shape function to interpolate the displacement field.
function plot_vibration_mode(SdivStruct, FemSol, referenceLength)
% PLOT_VIBRATION_MODE  Get an overview of the vibration modes.
%
% Arguments:
%	SdivStruct      (struct) -- Subdivised structure.
%	FemSol          (struct) -- Solution of the FEM simulation.
%	referenceLength (double) -- Reference length to scale the mode.

figure("WindowStyle", "docked");

% Scale factor.
filterTranslationDof = reshape((1:3)' + Node.nDof * ((1:SdivStruct.nNode)-1), 1, []);
maximumDeformation = @(iMode) max(abs(FemSol.mode(filterTranslationDof, iMode)));
percentageDeformation = 15;  % This gives a readable deformation.
computeScale = @(iMode) referenceLength/maximumDeformation(iMode) * percentageDeformation*1e-2;

% Plotting grid dimensions.
maxCols = 4;
if FemSol.nMode < maxCols
	nCols = FemSol.nMode;
else
	nCols = maxCols;
end
nRows = ceil(FemSol.nMode/nCols);

for i = 1:FemSol.nMode
	subplot(nRows, nCols, i);
	hold on;

	scale = computeScale(i);
	for elem = SdivStruct.elemList
		elem{:}.plotElem();  % pass `'Color', [0, 0, 0, 0.2]` for better clarity.
		x = [elem{:}.n1.pos(1) + scale * FemSol.mode(elem{:}.n1.dof(1), i),...
			 elem{:}.n2.pos(1) + scale * FemSol.mode(elem{:}.n2.dof(1), i)];
		y = [elem{:}.n1.pos(2) + scale * FemSol.mode(elem{:}.n1.dof(2), i), ...
			 elem{:}.n2.pos(2) + scale * FemSol.mode(elem{:}.n2.dof(2), i)];
		z = [elem{:}.n1.pos(3) + scale * FemSol.mode(elem{:}.n1.dof(3), i), ...
			 elem{:}.n2.pos(3) + scale * FemSol.mode(elem{:}.n2.dof(3), i)];
		plot3(x, y, z, Color=[0.9290 0.6940 0.1250], LineWidth=2);
	end

	xlabel("X/m");
	ylabel("Y/m");
	zlabel("Z/m");
	title(['f = ', num2str(FemSol.frequencyHertz(i)), ' Hz']);
	axis equal;
	grid;
	view(-35, 50);
	hold off;
end
end

%% 5. Total mass and sanity checks

function varargout = check_rbm(M_free, nNode, mass)
% CHECK_RBM  Perform sanity check based on rbm movement.
%
% This function uses a translation along the X-axis as a rigid body
% mode, to check that the total mass of the structure is the same as the
% one calculated with this RBM.
%
% Arguments:
%	M_free (nDof x nDof double) -- Global free mass matrix.
%	nNode  (int)                -- Number of structural nodes.
%	mass   (double)             -- Mass of the entire structure [kg].
% Return:
%	massFromRbm (double, optional) -- Mass calculated from RBM [kg].

% Rigid translation of 1m along the X-axis.
rbm = repmat([1, 0, 0, 0, 0, 0]', nNode, 1);

% Mass calculated from this translation.
massFromRbm = rbm' * M_free * rbm;

% Sanity checks.
allclose(mass, massFromRbm);

% Return the calculated mass, if wanted.
varargout = cell(nargout, 1);
for k = 1:nargout
	varargout{k} = massFromRbm;
end
end
