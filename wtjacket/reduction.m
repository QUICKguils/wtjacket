function ReductionSol = reduction(RunArg, Stm, SdivStruct, AlgSys)
% REDUCTION  Reduced models of the wt jacket.
%
% Arguments:
%	RunArg       (struct)     -- Code execution parameters, with fields:
%	  nMode      (int)        -- Number of first mode computed.
%	  tSet       (1xN double) -- Time sample used for time evolutions.
%	  nodeLabels (1xN double) -- Label list of nodes to inspect.
%	  m          (int)        -- Number of modes used in reductions.
%	  opts       (1xN char)   -- Output options.
%	    'p' -> Enable [P]lots creation.
%	Stm        (struct) -- Project statement data.
%	SdivStruct (struct) -- Subdivised structure.
%	AlgSys     (struct) -- Parameters of the discrete algebraic system.
% Returns:
%	ReductionSol (struct) -- Results of the reduction methods, with fields:
%	  ReducedSdivStruct (struct) -- Subdivised structure.
%	  CBReducedAlgSys   (struct) -- Parameters of the Craig-Brampton reduced discrete algebraic system.
%	  CBReducedFemSol   (struct) -- Solution of the reduced FEM simulation.
%	  ReducedNewmarkSol (struct) -- Solution of the reduced transient problem using newmark.

% Unpack relevant execution parameters.
LocalRunArg = {RunArg.nMode, RunArg.tSet, RunArg.nodeLabels, RunArg.m, RunArg.opts};
[nMode, tSet, nodeLabels, m, opts] = LocalRunArg{:};

% Remaining DOF per nodes.
%          x    y    z    \th_x \th_y \th_z
dofMask = [true true true false false true];

% Sort DOFs
nClampedDOFs = 24; % TODO improve this
[remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, nodeLabels, dofMask);

% 1. Guyans-Irons method - STATIC CONDENSATION
[~, GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol] = R_GuyanIronsReduction(SdivStruct, AlgSys, nMode, remainingDOFs, condensedDOFs);

% 2 - Craig-Brampton method
[CB_R, CBReducedSdivStruct, CBReducedAlgSys, CBReducedFemSol] = R_CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, remainingDOFs, condensedDOFs);

% Setup reduced load
[DiscreteLoad, FullLoad] = reducedLoad(CBReducedSdivStruct, SdivStruct, AlgSys, tSet, CB_R, remainingDOFs, condensedDOFs);

% Integrate
ReducedNewmarkSol = newmark_integrate_reduced_system(Stm, CBReducedSdivStruct, CBReducedAlgSys, tSet, DiscreteLoad, FullLoad, nMode, opts);

% 4. Build the return data structure

ReductionSol.GIReducedSdivStruct = GIReducedSdivStruct;
ReductionSol.GIReducedAlgSys = GIReducedAlgSys;
ReductionSol.GIReducedFemSol = GIReducedFemSol;
ReductionSol.CBReducedAlgSys = CBReducedAlgSys;
ReductionSol.CBReducedFemSol = CBReducedFemSol;
ReductionSol.ReducedNewmarkSol = ReducedNewmarkSol;
end

function [DiscreteLoad, FullLoad] = reducedLoad(CBReducedSdivStruct, SdivStruct, AlgSys, tSet, R, remainingDOFs, condensedDOFs)
% Choose the load to study.
lookupLoadLabel = 1;
FullLoad = SdivStruct.loadList{lookupLoadLabel};
FullDiscreteLoad = FullLoad.set_discrete_load(AlgSys.nDofFree, tSet);
s = FullDiscreteLoad.sample;

% Ignore clamped dofs
s = s(~AlgSys.cstrMask, :);

% Reorder DOFs
S = zeros(size(s));
for t =1:length(tSet)
	S(:, t) = [s(remainingDOFs, t); s(condensedDOFs, t)];
end

% Reduce discretized applied load
for t = 1:length(tSet)
	DiscreteLoad.sample(:, t) = R' * S(:, t); % using R matrix from CB
end
DiscreteLoad.node = CBReducedSdivStruct.nodeList{1, 1};
end

function ReducedNewmarkSol = newmark_integrate_reduced_system(Stm, ReducedSdivStruct, ReducedAlgSys, tSet, DiscreteLoad, ThisLoad, nMode, opts)
% NEWMARK INTEGRATE REDUCED SYSTEM

% Time the time structure that Newmak need.
TimeParams.sample = tSet;
TimeParams.steps  = [diff(tSet), tSet(end)-tSet(end-1)];
TimeParams.numel  = numel(tSet);
TimeParams.initialConditions = Stm.INITIAL_CONDITIONS;

% Solve with newmark.
ReducedNewmarkSol = newmark(ReducedAlgSys, TimeParams, DiscreteLoad.sample);
ReducedNewmarkSol.TimeParams = TimeParams;
% Plot displacements at hightlighted nodes.
if contains(opts, 'p')
	lookupNodeLabels = [1, 2];

	plot_displacement(ReducedNewmarkSol, tSet, lookupNodeLabels, ThisLoad.direction, ReducedSdivStruct.nodeList, nMode);
end

	function plot_displacement(TransientSol, tSet, lookupNodeLabels, loadDirection, nodeList, nMode)
		% PLOT_DISPLACEMENT  Plot the time evolution of the displacements.

		allclose(norm(loadDirection), 1);

		figure("WindowStyle", "docked");

		nNode = numel(lookupNodeLabels);
		for i = 1:nNode
			qX = TransientSol.q(nodeList{lookupNodeLabels(i)}.dof(1), :) * loadDirection(1);
			qY = TransientSol.q(nodeList{lookupNodeLabels(i)}.dof(2), :) * loadDirection(2);
			qZ = TransientSol.q(nodeList{lookupNodeLabels(i)}.dof(3), :) * loadDirection(3);

			qDir = qX + qY + qZ;

			subplot(nNode, 1, i);
			plot(tSet, qDir);
			xlabel("Time (s)");
			ylabel("Displacement (dir: [" + num2str(loadDirection, '%.3f  ') + "])");
			title('Transient response', ...
				['(node: ', num2str(lookupNodeLabels(i)), ...
				', method: ', TransientSol.name, ...
				', order: ', num2str(nMode), ')']);
			grid;
		end
	end

end

function [remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, nodeLabels, dofMask)
% SORT DOFS Sort structural dofs into remainingDOFs, condensedDOFs vectors
% Arguments:
%	SdivStruct          (struct) -- Subdivised structure.
%	AlgSys              (struct) -- Parameters of the discrete algebraic system.
%	highlightedNodes    (nNodes x Node) -- Array of Node Index used for sorting.
%	dofMask             (6 x boolean) -- True for kept dof at nodes (x, y, z, theta_x, theta_y, theta_z)
% Returns:
%   remainingDOFs       (nb x int) -- Array of selected remaining dofs of the structure
%   condensedDOFs       (nDOF-nb x int) -- Array of condensed excluded dofs of the structure

remainingDOFs = [];
for k=1:length(nodeLabels)
	remainingDOFs = [remainingDOFs SdivStruct.nodeList{1, nodeLabels(k)}.dof(dofMask)];
end

% since clamped dofs
remainingDOFs = remainingDOFs - 24;

condensedDOFs = zeros(1, AlgSys.nDof-length(remainingDOFs));
k = 1;
l = 1;
while k<=AlgSys.nDof
	if not(ismember(k, remainingDOFs))
		condensedDOFs(l) = k;
		l=l+1;
	end
	k=k+1;
end
end

function ReducedFemSol = solve_eigenvalue_problem_mass_normalized(AlgSys, nMode)
% SOLVE EIGENVALUE PROBLEM MASS NORMALIZED
% Arguments:
%   Algsys
%   nMode
% Returns:
%   ReducedFemSol

[Modes, ReducedFemSol.Omega] = eigs(AlgSys.K, AlgSys.M, nMode, 'sm');
% Mass-normalize modes per definition
for i = 1:nMode
	MassNormalizedModes(:, i) = Modes(:, i) / sqrt(Modes(:, i)' * AlgSys.M * Modes(:, i));
end
ReducedFemSol.mode = MassNormalizedModes;
ReducedFemSol.nMode = nMode;
ReducedFemSol.frequencyRad = sqrt(diag(ReducedFemSol.Omega));
ReducedFemSol.frequencyHertz = ReducedFemSol.frequencyRad/(2*pi);
end

function [reorderedA, Arr, Arc, Acr, Acc] = submatrixRecomposition(A, r, c)
% SUBMATRIX RECOMPOSITION Decomposes A into Arr, Arc, Acr, Acc, given r, c
% A = [Arr Arc; Acr Acc]
% Arguments:
%   A   matrix
%   r   index for r part of decomposition
%   c   index for c part of decomposition
% Returns:
%   reorderedA
%   Arr
%   Arc
%   Acr
%   Acc

Arr = A(r, r);
Arc = A(r, c);
Acr = A(c, r);
Acc = A(c, c);

% Reordering
reorderedA = [Arr Arc; Acr Acc];
end

function ReducedAlgSys = ReduceAlgSys(AlgSys, R)
% REDUCE ALGSYS Reduces AlgSys.K/M/C with R

ReducedAlgSys.K = R' * AlgSys.K * R;
ReducedAlgSys.M = R' * AlgSys.M * R;
ReducedAlgSys.C = R' * AlgSys.C * R;

ReducedAlgSys.M_free = ReducedAlgSys.M;
ReducedAlgSys.K_free = ReducedAlgSys.K;
end

% 1 - Guyan Irons Reduction
function [R, ReducedSdivStruct, ReducedAlgSys, ReducedFemSol] = R_GuyanIronsReduction(SdivStruct, AlgSys, nMode, remainingDOFs, condensedDOFs)
ReducedSdivStruct.nodeList{1, 1} = SdivStruct.nodeList{1, 18};
ReducedSdivStruct.nodeList{1, 2} = SdivStruct.nodeList{1, 22};

ReducedSdivStruct.nodeList{1, 1}.dof = [1:4 0 0];
ReducedSdivStruct.nodeList{1, 2}.dof = [5:8 0 0];

K = AlgSys.K;
M = AlgSys.M;
C = AlgSys.C;

% Reordering matrices for reduction.
[AlgSys.K, ~, ~, KCR, KCC] = submatrixRecomposition(K, remainingDOFs, condensedDOFs);
[AlgSys.M, ~, ~,   ~,   ~] = submatrixRecomposition(M, remainingDOFs, condensedDOFs);
[AlgSys.C, ~, ~,   ~,   ~] = submatrixRecomposition(C, remainingDOFs, condensedDOFs);

% Reduction matrix R computation.
R = [eye(length(remainingDOFs)); -inv(KCC) * KCR];

% ReducedAlgSys computation.
ReducedAlgSys = ReduceAlgSys(AlgSys, R);

% Ensure remaining dofs are free
ReducedAlgSys.cstrMask = false(1, length(remainingDOFs));

ReducedAlgSys.nDof = length(remainingDOFs);
ReducedAlgSys.nDofFree = ReducedAlgSys.nDof;

% Solve the eigenvalue problem.
ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end

% 2 - Craig-Brampton Reduction
function  [R, ReducedSdivStruct, ReducedAlgSys, ReducedFemSol] = R_CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, B, I)
% R CRAIG BRAMPTON REDUCTION Reduce the AlgSys using the R matrix reduction, using m submodes


ReducedSdivStruct.nodeList{1, 1} = SdivStruct.nodeList{1, 18};
ReducedSdivStruct.nodeList{1, 2} = SdivStruct.nodeList{1, 22};

ReducedSdivStruct.nodeList{1, 1}.dof = [1:4 0 0];
ReducedSdivStruct.nodeList{1, 2}.dof = [5:8 0 0];

K = AlgSys.K;
M = AlgSys.M;
C = AlgSys.C;

% Reordering matrices
[AlgSys.K, ~, ~, KIB, KII] = submatrixRecomposition(K, B, I);
[AlgSys.M, ~, ~,   ~, MII] = submatrixRecomposition(M, B, I);
[AlgSys.C, ~, ~,   ~,   ~] = submatrixRecomposition(C, B, I);

% Solving substructure fixed at the boundary - for m modes, mass-normalized
SubAlgSys.M = MII;
SubAlgSys.K = KII;
SubFemSol = solve_eigenvalue_problem_mass_normalized(SubAlgSys, m);

% R Matrix computation
nb = length(B);
R = [eye(nb, nb+m); -inv(KII)*KIB SubFemSol.mode];

% ReducedAlgSys computation.
ReducedAlgSys = ReduceAlgSys(AlgSys, R);

ReducedAlgSys.cstrMask = false(1, length(B)+m);

ReducedAlgSys.nDof = length(B)+m;
ReducedAlgSys.nDofFree = ReducedAlgSys.nDof;

% Solve the reduced eigenvalue problem.
ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end
