function [GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol, CBReducedAlgSys, CBReducedFemSol, ReducedNewmarkSol] = reduction(Cst, SdivStruct, AlgSys, nMode, m, opts)
% REDUCTION Study of the reduced models of the wt jacket, from full model.
%
% Arguments:
%	Cst         (struct) -- Constant project quantities.
%	SdivStruct  (struct) -- Subdivised structure.
%	AlgSys      (struct) -- Parameters of the discrete algebraic system.
%	nMode       (int)    -- Number of computed modes.
%	opts  (1xN char) -- Options.
%	  ''  -> No options.
%	  'p' -> Enable plots creation.
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

% Sort DOFs
nClampedDOFs = 24; % TODO improve this
[remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, highlightedNodes, dofMask, nClampedDOFs);

% 1 - Guyans-Irons Method - STATIC CONDENSATION.
[GI_R, GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol] = R_GuyanIronsReduction(SdivStruct, AlgSys, nMode, remainingDOFs, condensedDOFs);

% 2 - Craig-Brampton Method.
[CB_R, CBReducedSdivStruct, CBReducedAlgSys, CBReducedFemSol] = R_CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, remainingDOFs, condensedDOFs);

% 3 - Time integration.
timeSample = 0:0.01:10;
TimeParams = set_time_parameters(timeSample, Cst.INITIAL_CONDITIONS);

% Setup reduced load
[DiscreteLoad, FullLoad] = reducedLoad(CBReducedSdivStruct, SdivStruct, AlgSys, TimeParams, CB_R, remainingDOFs, condensedDOFs);

% Integrate
ReducedNewmarkSol = newmark_integrate_reduced_system(Cst, CBReducedSdivStruct, CBReducedAlgSys, TimeParams, DiscreteLoad, FullLoad, nMode, opts);
end

function [DiscreteLoad, FullLoad] = reducedLoad(CBReducedSdivStruct, SdivStruct, AlgSys, TimeParams, R, remainingDOFs, condensedDOFs)
% Choose the load to study.
lookupLoadLabel = 1;
FullLoad = SdivStruct.loadList{lookupLoadLabel};

FullDiscreteLoad = FullLoad.set_discrete_load(AlgSys.nDof_free, TimeParams.sample);
s = FullDiscreteLoad.sample;

% Ignore clamped dofs
s = s(~AlgSys.cstrMask, :);

% Reorder DOFs
S = zeros(size(s));
for t =1:length(TimeParams.sample)
    S(:, t) = [s(remainingDOFs, t); s(condensedDOFs, t)];
end

% Reduce discretized applied load
for t = 1:length(TimeParams.sample)
    DiscreteLoad.sample(:, t) = R' * S(:, t); % using R matrix from CB
end
DiscreteLoad.node = CBReducedSdivStruct.nodeList{1, 1};
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

function [remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, highlightedNodes, dofMask, nClampedDOFs)
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
for k=1:length(highlightedNodes)
    remainingDOFs = [remainingDOFs SdivStruct.nodeList{1, highlightedNodes(k)}.dof(dofMask)];
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
ReducedAlgSys.nDof_free = ReducedAlgSys.nDof;

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
ReducedAlgSys.nDof_free = ReducedAlgSys.nDof;

% Solve the reduced eigenvalue problem.
ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end