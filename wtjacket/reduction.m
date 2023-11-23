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
    [remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, highlightedNodes, dofMask);

    % 1 - Guyans-Irons Method - STATIC CONDENSATION.
    [GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol] = R_GuyanIronsReduction(SdivStruct, AlgSys, nMode, remainingDOFs, condensedDOFs, highlightedNodes);

    % 2 - Craig-Brampton Method.
    [CBReducedSdivStruct, CBReducedAlgSys, CBReducedFemSol] = R_CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, remainingDOFs, condensedDOFs, highlightedNodes);
    CBReducedFemSol.frequencyHertz
    % 3 - Time integration.
    timeSample = 0:0.1:10;
    TimeParams = set_time_parameters(timeSample, Cst.INITIAL_CONDITIONS);

    ReducedNewmarkSol = newmark_integrate_reduced_system(Cst, CBReducedSdivStruct, CBReducedAlgSys, CBReducedFemSol, TimeParams, nMode, opts);
end

function ReducedNewmarkSol = newmark_integrate_reduced_system(Cst, ReducedSdivStruct, ReducedAlgSys, ReducedFemSol, TimeParams, nMode, opts)
    % Applied load.
    loadAmplitude = Cst.TAIL_MASS * Cst.TAIL_SPEED * Cst.MOMENTUM_TRANSFER / Cst.IMPACT_DURATION;
    loadDirection = [cosd(Cst.FORCE_DIRECTION), -sind(Cst.FORCE_DIRECTION), 0];

    ThisLoad = HarmonicLoad(ReducedSdivStruct.nodeList{1, 1}, loadDirection, loadAmplitude, Cst.LOAD_FREQUENCY_HERTZ);
    
    % Create the time-discretized load.
    DiscreteLoad = ThisLoad.set_discrete_load(ReducedAlgSys.nDof, TimeParams.sample);

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
    for k=1:length(highlightedNodes)
        remainingDOFs = [remainingDOFs SdivStruct.nodeList{1, highlightedNodes(k)}.dof(dofMask)];
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

% 1 - Guyan Irons Reduction

function [ReducedSdivStruct, ReducedAlgSys, ReducedFemSol] = R_GuyanIronsReduction(SdivStruct, AlgSys, nMode, remainingDOFs, condensedDOFs, highlightedNodes)
    % TODO refactor as suppressed collumns ?
    ReducedSdivStruct.nodeList{1, 1} = SdivStruct.nodeList{1, 18};
    ReducedSdivStruct.nodeList{1, 2} = SdivStruct.nodeList{1, 22};
    
    ReducedSdivStruct.nodeList{1, 1}.dof = [1:4 0 0];
    ReducedSdivStruct.nodeList{1, 2}.dof = [5:8 0 0];

    % K Submatrices
    KCC = AlgSys.K(condensedDOFs, condensedDOFs);
    KCR = AlgSys.K(condensedDOFs, remainingDOFs);
    KRR = AlgSys.K(remainingDOFs, remainingDOFs);
    KRC = AlgSys.K(remainingDOFs, condensedDOFs);

    % M Submatrices
    MCC = AlgSys.M(condensedDOFs, condensedDOFs);
    MCR = AlgSys.M(condensedDOFs, remainingDOFs);
    MRR = AlgSys.M(remainingDOFs, remainingDOFs);
    MRC = AlgSys.M(remainingDOFs, condensedDOFs);

    % C Submatrices
    CCC = AlgSys.C(condensedDOFs, condensedDOFs);
    CCR = AlgSys.C(condensedDOFs, remainingDOFs);
    CRR = AlgSys.C(remainingDOFs, remainingDOFs);
    CRC = AlgSys.C(remainingDOFs, condensedDOFs);
    
    % Reordering
    ReducedAlgSys.M = [MRR MRC; MCR MCC];
    ReducedAlgSys.K = [KRR KRC; KCR KCC];
    ReducedAlgSys.C = [CRR CRC; CCR CCC];
    
    % R Matrix computation
    R = [eye(max(size(remainingDOFs))); -inv(KCC) * KCR];
    
    % ReducedAlgSys computation.
    ReducedAlgSys.M = R' * ReducedAlgSys.M * R;
    ReducedAlgSys.K = R' * ReducedAlgSys.K * R;
    ReducedAlgSys.C = R' * ReducedAlgSys.C * R;
    
    ReducedAlgSys.M_free = ReducedAlgSys.M;
    ReducedAlgSys.K_free = ReducedAlgSys.K;

    ReducedAlgSys.cstrMask = false(1, length(remainingDOFs));

    ReducedAlgSys.nDof = length(remainingDOFs);
    ReducedAlgSys.nDof_free = ReducedAlgSys.nDof;

    % Solve the eigenvalue problem.
    ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end

% 2 - Craig-Brampton Reduction

function  [ReducedSdivStruct, ReducedAlgSys, ReducedFemSol] = R_CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, remainingDOFs, condensedDOFs, highlightedNodes)
    ReducedSdivStruct.nodeList{1, 1} = SdivStruct.nodeList{1, 18};
    ReducedSdivStruct.nodeList{1, 2} = SdivStruct.nodeList{1, 22};
    
    ReducedSdivStruct.nodeList{1, 1}.dof = [1:4 0 0];
    ReducedSdivStruct.nodeList{1, 2}.dof = [5:8 0 0];
    
    % K Submatrices
    KII = AlgSys.K(condensedDOFs, condensedDOFs);
    KIB = AlgSys.K(condensedDOFs, remainingDOFs);
    KBB = AlgSys.K(remainingDOFs, remainingDOFs);
    KBI = AlgSys.K(remainingDOFs, condensedDOFs);    
    
    % M Submatrices
    MII = AlgSys.M(condensedDOFs, condensedDOFs);
    MIB = AlgSys.M(condensedDOFs, remainingDOFs);
    MBB = AlgSys.M(remainingDOFs, remainingDOFs);
    MBI = AlgSys.M(remainingDOFs, condensedDOFs);

    % C Submatrices
    CII = AlgSys.C(condensedDOFs, condensedDOFs);
    CIB = AlgSys.C(condensedDOFs, remainingDOFs);
    CBB = AlgSys.C(remainingDOFs, remainingDOFs);
    CBI = AlgSys.C(remainingDOFs, condensedDOFs);

    % Reordering
    ReducedAlgSys.M = [MBB MBI; MIB MII];
    ReducedAlgSys.K = [KBB KBI; KIB KII];
    ReducedAlgSys.C = [CBB CBI; CIB CII];

    % Solving substructure fixed at the boundary - for m modes, mass-normalized
    SubAlgSys.M = MII;
    SubAlgSys.K = KII;
    SubFemSol = solve_eigenvalue_problem_mass_normalized(SubAlgSys, m);
    
    % R Matrix computation
    nb = length(remainingDOFs);
    R = [eye(nb) zeros(nb, m); -inv(KII) * KIB SubFemSol.mode];
    
    % ReducedAlgSys computation.
    ReducedAlgSys.M = R' * ReducedAlgSys.M * R;
    ReducedAlgSys.K = R' * ReducedAlgSys.K * R;
    ReducedAlgSys.C = R' * ReducedAlgSys.C * R;

    ReducedAlgSys.M_free = ReducedAlgSys.M;
    ReducedAlgSys.K_free = ReducedAlgSys.K;

    ReducedAlgSys.cstrMask = false(1, length(remainingDOFs)+m);

    ReducedAlgSys.nDof = length(remainingDOFs) + m;
    ReducedAlgSys.nDof_free = ReducedAlgSys.nDof;

    % Solve the reduced eigenvalue problem.
    ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end