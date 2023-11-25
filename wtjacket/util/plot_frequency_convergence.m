function plot_frequency_convergence(Solution)
    % PLOT_FREQUENCY_CONVERGENCE  Generate convergence graphs for the frequencies.
    %
    % Arguments:
    %	Solution (struct) -- Set of solution, with fields:
    %	  nMode   (int)        -- Number of first mode computed.
    %	  sdivSet (1xN int)    -- Set of number of subdivisions.
    %	  freqSet (1xN double) -- Set of associated frequencies.
    %	  name    (char)       -- Name of the solution.
    
    % Build the solution name.
    if Solution.name ~= ""
        Solution.name = " (" + Solution.name + ")";
    end
    
    % Instantiate a figure object.
    figure("WindowStyle", "docked");
    
    % Plot the frequencies.
    subplot(1, 2, 1);
    plot(Solution.sdivSet, Solution.freqSet);
    title("Frequencies convergence" + Solution.name);
    xlabel("Number of sub-elements");
    ylabel("Natural frequency (Hz)");
    grid;
    
    % Plot the residuals.
    subplot(1, 2, 2);
    % NOTE:
    % this way of defining residuals is taken from the fluid mechanics
    % course, section 3.3.
    residuals = abs(diff(Solution.freqSet)) ./ abs(Solution.freqSet(2, :) - Solution.freqSet(1, :));
    length(residuals)
    semilogy(Solution.sdivSet(1:end-1), residuals);
    title("Residuals convergence" + Solution.name);
    xlabel("Number of sub-elements");
    ylabel("Residual");
    grid;
end
