function reduction_analysis(Cst, SdivStruct, AlgSys, nMode, opts)
    % REDUCTION_ANALYSIS  Analyze the results of the reduction part.
    % TODO - freq convergence plot
    m = 1
    [ReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol, CBReducedAlgSys, CBReducedFemSol, ReducedNewmarkSol] = reduction(Cst, SdivStruct, AlgSys, nMode, m, opts);
    FrequencyHertzArray = [GIReducedFemSol.frequencyHertz CBReducedFemSol.frequencyHertz];
    mModesInCB = 5:5:50
    for m=1:length(mModesInCB)
        [ReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol, CBReducedAlgSys, CBReducedFemSol, ReducedNewmarkSol] = reduction(Cst, SdivStruct, AlgSys, nMode, mModesInCB(m), opts);
        FrequencyHertzArray = [FrequencyHertzArray CBReducedFemSol.frequencyHertz];
    end
    Sol.nMode = nMode
    Sol.sdivSet = [0 1 mModesInCB]
    Sol.freqSet = FrequencyHertzArray
    Sol.name = "Craig Brampton Reduced Model"
    %	  sdivSet (1xN int)    -- Set of number of subdivisions.
    %	  freqSet (1xN double) -- Set of associated frequencies.
    %	  name    (char)       -- Name of the solution.

    plot_frequency_convergence(Sol)
end