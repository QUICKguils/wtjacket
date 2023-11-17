function [ReducedSdivStruct, CBReducedAlgSys, CBReducedFemSol, ReducedNewmarkSol] = reduction(Cst, SdivStruct, AlgSys, nMode, opts)
    % TODO documentation
    % REDUCTION  

    % Nodes highlighted in the structure
    highlightedNodes = [18, 22];

    % Remaining DOF per nodes
    % x y z \theta_z
    dofMask = [true true true false false true]; 
    
    % Sort DOFs
    [remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, highlightedNodes, dofMask);

    % Guyans-Irons Method - STATIC CONDENSATION
    [GIReducedAlgSys, GIReducedFemSol] = GuyanIronsReduction(AlgSys, nMode, remainingDOFs, condensedDOFs, highlightedNodes);

    % Craig-Brampton Method 
    % TODO refactor this part into analysis/reduction_analysis
    m = 10;
    [ReducedSdivStruct, CBReducedAlgSys, CBReducedFemSol] = CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, remainingDOFs, condensedDOFs, highlightedNodes);
    % TODO - freq convergence plot

    % 3
    timeSample = 0:0.01:10;
    TimeParams = set_time_parameters(timeSample, Cst.INITIAL_CONDITIONS);

    % C matrix
    eps = [Cst.DAMPING_RATIO, Cst.DAMPING_RATIO];
    [CBReducedAlgSys.C, eps] = set_damping_parameters(eps, CBReducedFemSol.frequencyRad, CBReducedAlgSys);

    % Choose the load to study.
    loadAmplitude = Cst.TAIL_MASS * Cst.TAIL_SPEED * Cst.MOMENTUM_TRANSFER / Cst.IMPACT_DURATION;
    loadDirection = [cosd(Cst.FORCE_DIRECTION), -sind(Cst.FORCE_DIRECTION), 0];
    
    SdivStruct.nodeList{1, 18}.dof = [1:4 0 0];
    SdivStruct.nodeList{1, 22}.dof = [5:8 0 0];
    ThisLoad = HarmonicLoad(SdivStruct.nodeList{1, 18}, loadDirection, loadAmplitude, Cst.LOAD_FREQUENCY_HERTZ);
    % Create the time-discretized load.
    DiscreteLoad = ThisLoad.set_discrete_load(CBReducedAlgSys.nDof, TimeParams.sample);

    % Solve with newmark
    ReducedNewmarkSol = newmark(CBReducedAlgSys, TimeParams, DiscreteLoad.sample);
    ReducedNewmarkSol.q(:,100);
    % Plot
    if contains(opts, 'p')
        lookupNodeLabels = highlightedNodes;
        plot_displacement(ReducedNewmarkSol, TimeParams.sample, lookupNodeLabels, ThisLoad.direction, ReducedSdivStruct.nodeList, nMode);
    end
end

function [remainingDOFs, condensedDOFs] = sort_DOFS(SdivStruct, AlgSys, highlightedNodes, dofMask)
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

function [ReducedAlgSys, ReducedFemSol] = GuyanIronsReduction(AlgSys, nMode, remainingDOFs, condensedDOFs, highlightedNodes)
    % TODO refactor as suppressed collumns ? 

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

    % ReducedAlgSys computation.
    ReducedAlgSys.M = MRR - MRC / KCC * KCR - KRC / KCC * MCR + KRC / KCC * MCC / KCC * KCR;
    ReducedAlgSys.K = KRR - KRC / KCC * KCR;
    
    ReducedAlgSys.M_free = ReducedAlgSys.M;
    ReducedAlgSys.K_free = ReducedAlgSys.K;

    ReducedAlgSys.nDof = length(remainingDOFs);
    % Solve the eigenvalue problem.
    ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end

% 2 - Craig-Brampton Reduction

function  [ReducedSdivStruct, ReducedAlgSys, ReducedFemSol] = CraigBramptonReduction(SdivStruct, AlgSys, nMode, m, remainingDOFs, condensedDOFs, highlightedNodes)
    ReducedSdivStruct.nodeList{1, 18}.dof = [1:4];
    ReducedSdivStruct.nodeList{1, 22}.dof = [5:8];
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

    SubAlgSys.M = MII;
    SubAlgSys.K = KII;
    % Solving substructure fixed at the boundary - for m modes, mass-normalized
    SubFemSol = solve_eigenvalue_problem_mass_normalized(SubAlgSys, m);
    
    % Building reduced K.
    KBBtilde = KBB - KBI / KII * KIB;
    K_reduced = [KBBtilde zeros(length(remainingDOFs), m); zeros(m, length(remainingDOFs)) SubFemSol.Omega];

    % Building reduced M.
    MBBtilde = MBB - (MBI/KII) * MIB + (KBI/KII) * (MII/KII) * KIB;
    XI = SubFemSol.mode;
    MBtilde = XI' * (MIB - (MII / KII) * KIB);
    M_reduced = [MBBtilde MBtilde'; MBtilde eye(m)];

    % ReducedAlgSys computation.
    ReducedAlgSys.M = M_reduced;
    ReducedAlgSys.K = K_reduced;

    ReducedAlgSys.M_free = ReducedAlgSys.M;
    ReducedAlgSys.K_free = ReducedAlgSys.K;

    ReducedAlgSys.cstrMask = repmat(false, 1, length(remainingDOFs) + m);

    ReducedAlgSys.nDof = length(remainingDOFs) + m;
    ReducedAlgSys.nDof_free = ReducedAlgSys.nDof;

    % Solve the reduced eigenvalue problem.
    ReducedFemSol = solve_eigenvalue_problem_mass_normalized(ReducedAlgSys, nMode);
end