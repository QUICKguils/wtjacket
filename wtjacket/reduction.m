function [GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol, CBReducedAlgSys, CBReducedFemSol, ReducedNewmarkSol] = reduction(Cst, SdivStruct, AlgSys, nMode, m, opts)
% REDUCTION Study of the reduced models of the wt jacket, from full model.
%
% Arguments:
%	Cst         (struct)  -- Constant project quantities.
%	SdivStruct  (struct)  -- Subdivised structure.
%	AlgSys      (struct)  -- Parameters of the discrete algebraic system.
%	nMode       (int)     -- Number of computed modes.
%	opts       (1xN char) -- Options.
%	  'p' -> Enable [P]lots creation.
% Returns:
%	ReducedSdivStruct  (struct) -- Subdivised structure.
%	CBReducedAlgSys      (struct) -- Parameters of the Craig-Brampton reduced discrete algebraic system.
%	CBReducedFemSol      (struct) -- Solution of the reduced FEM simulation.
%	ReducedNewmarkSol (struct) -- Solution of the reduced transient problem using newmark.

% Nodes highlighted in the structure.
highlightedNodes = [18, 22];

% Remaining DOF per nodes.
%          x    y    z    \th_x \th_y \th_z
dofMask = [true true true false false true];

% Sort DOFs.
[remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, inspectNodeLabels, dofMask);

% 1 - Guyans-Irons Method - STATIC CONDENSATION.
[GI_R, GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol] = R_GuyanIronsReduction(SdivStruct, AlgSys, nMode, remainingDOFs, condensedDOFs, highlightedNodes);
% 2 - Craig-Brampton Method.
[CB_R, CBReducedSdivStruct, CBReducedAlgSys, CBReducedFemSol] = R_CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, remainingDOFs, condensedDOFs, highlightedNodes);

% 3 - Time integration.
timeSample = 0:0.01:10;
TimeParams = set_time_parameters(timeSample, Cst.INITIAL_CONDITIONS);

% Applied load.
loadAmplitude = Cst.TAIL_MASS * Cst.TAIL_SPEED * Cst.MOMENTUM_TRANSFER / Cst.IMPACT_DURATION;
loadDirection = [cosd(Cst.FORCE_DIRECTION), -sind(Cst.FORCE_DIRECTION), 0];

ThisLoad = HarmonicLoad(CBReducedSdivStruct.nodeList{1, 1}, loadDirection, loadAmplitude, Cst.LOAD_FREQUENCY_HERTZ);

% Create the time-discretized load.
DiscreteLoad = ThisLoad.set_discrete_load(CBReducedAlgSys.nDof_free, TimeParams.sample);

% Choose the load to study.
lookupLoadLabel = 1;
FullLoad = SdivStruct.loadList{lookupLoadLabel};
FullDiscreteLoad = FullLoad.set_discrete_load(AlgSys.nDof_free, TimeParams.sample);
s = FullDiscreteLoad.sample;
s = s(~AlgSys.cstrMask, :);

S = zeros(size(s));
for t =1:length(timeSample)
    S(:, t) = [s(remainingDOFs, t); s(condensedDOFs, t)];
end
for t = 1:length(timeSample)
    DiscreteLoad.sample(:, t) = CB_R' * S(:, t);
end
DiscreteLoad.node = CBReducedSdivStruct.nodeList{1, 1};

% Reduce applied load

ReducedNewmarkSol = newmark_integrate_reduced_system(Cst, CBReducedSdivStruct, CBReducedAlgSys, TimeParams, DiscreteLoad, ThisLoad, nMode, opts);
end

function ReducedNewmarkSol = newmark_integrate_reduced_system(Cst, ReducedSdivStruct, ReducedAlgSys, TimeParams, DiscreteLoad, ThisLoad, nMode, opts)
% NEWMARK INTEGRATE REDUCED SYSTEM
% Solve with newmark.
ReducedNewmarkSol = newmark(ReducedAlgSys, TimeParams, DiscreteLoad.sample);

% Plot displacements at hightlighted nodes.
if contains(opts, 'p')
    lookupNodeLabels = [1, 2];
    plot_displacement(ReducedNewmarkSol, TimeParams.sample, lookupNodeLabels, ThisLoad.direction, ReducedSdivStruct.nodeList, nMode);
end
end

function [remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, highlightedNodes, dofMask)
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
for k=1:length(inspectNodeLabels)
	remainingDOFs = [remainingDOFs SdivStruct.nodeList{1, inspectNodeLabels(k)}.dof(dofMask)];
end

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
function [R, ReducedSdivStruct, ReducedAlgSys, ReducedFemSol] = R_GuyanIronsReduction(SdivStruct, AlgSys, nMode, remainingDOFs, condensedDOFs, highlightedNodes)
% TODO refactor as suppressed collumns ?
ReducedSdivStruct.nodeList{1, 1} = SdivStruct.nodeList{1, 18};
ReducedSdivStruct.nodeList{1, 2} = SdivStruct.nodeList{1, 22};

ReducedSdivStruct.nodeList{1, 1}.dof = [1:4 0 0];
ReducedSdivStruct.nodeList{1, 2}.dof = [5:8 0 0];

K = AlgSys.K;
M = AlgSys.M;
C = AlgSys.C;

% Reordering matrices for reduction.
[AlgSys.K, ~, ~, KCR, KCC] = submatrixRecomposition(AlgSys.K, remainingDOFs, condensedDOFs);
[AlgSys.M, ~, ~,   ~,   ~] = submatrixRecomposition(AlgSys.M, remainingDOFs, condensedDOFs);
[AlgSys.C, ~, ~,   ~,   ~] = submatrixRecomposition(AlgSys.C, remainingDOFs, condensedDOFs);

% Reduction matrix R computation.
R = [eye(length(remainingDOFs)); -inv(KCC) * KCR];

% ReducedAlgSys computation.
ReducedAlgSys = ReduceAlgSys(AlgSys, R);

% Ensure remaining dofs are free
ReducedAlgSys.cstrMask = false(1, length(remainingDOFs));

ReducedAlgSys.nDof = length(remainingDOFs);
ReducedAlgSys.nDof_free = ReducedAlgSys.nDof;

% Solve the eigenvalue problem.
ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end

% 2 - Craig-Brampton Reduction
function  [R, ReducedSdivStruct, ReducedAlgSys, ReducedFemSol] = R_CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, B, I, highlightedNodes)
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
ReducedAlgSys.nDof_free = ReducedAlgSys.nDof;

% Solve the reduced eigenvalue problem.
ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end
