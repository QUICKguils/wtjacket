%% Plot the displacements

function plot_displacement(TransientSol, timeSample, lookupNodeLabels, loadDirection, nodeList, nMode)
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
        plot(timeSample, qDir);
        xlabel("Time (s)");
        ylabel("Displacement (dir: [" + num2str(loadDirection, '%.3f  ') + "])");
        title('Transient response', ...
            ['(node: ', num2str(lookupNodeLabels(i)), ...
            ', method: ', TransientSol.name, ...
            ', order: ', num2str(nMode), ')']);
        grid;
    end
end
