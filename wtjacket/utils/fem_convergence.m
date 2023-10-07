function fem_convergence
	% CONVERGENCE  Analyze the convergence of the FEM simulation.

	% TODO:
	% - analyse the NX modes.
	% - implement other ways to analyze the convergence, e.g. MAC w.r.t. the
	%   NX modes ?

    % TODO mass convergence analysis

	% Directory where the present file lies.
	file_dir = fileparts(mfilename("fullpath"));
	% Results file.
	res_file = (fullfile(file_dir, "../../res/"));

	% Run the bare structure simulation once, if SOL data are not available.
	if ~isfile(fullfile(res_file, "modeling_sol.mat"))
		modeling(1, '');
	end
	SOL = load(fullfile(res_file, "modeling_sol.mat"));
    BS = load(fullfile(res_file, "bare_struct.mat"));
    

	% Desired sdiv sequence.
	sdiv_serie = 1:4;
    
	freqs_serie = zeros(numel(sdiv_serie), SOL.nbMode);
    modes_serie = zeros(numel(sdiv_serie), BS.nbDOF, SOL.nbMode);    
    masses_serie = zeros(numel(sdiv_serie), 2);
    
	% Gather the frequencies for the desired subdivisions.
	for i = 1:numel(sdiv_serie)
		modeling(sdiv_serie(i), '');
		SOL = load("res\modeling_sol.mat");
        BS = load(fullfile(res_file, "bare_struct.mat"));

		masses_serie(i, 1) = SOL.mass_rbm;
        masses_serie(i, 2) = SOL.mass_rbm/BS.mass;
        
        freqs_serie(i, :) = SOL.freqs;
        
        for d = 1:BS.nbDOF
            modes_serie(i, d, :) = SOL.modes(d,:);
		end
    end
    
    % Mass convergence
    figure("WindowStyle", "docked");
    subplot(1,2,1);
    plot(sdiv_serie, masses_serie(:, 1));
	title("RBM Total Mass convergence");
	xlabel("Number of sub-elements");
	ylabel("Total mass from RBM (kg)");
	grid;
    
    subplot(1,2,2);
    plot(sdiv_serie, masses_serie(:, 2));
	title("RBM Total Mass - Theoretical Mass Ratio convergence"); % TODO better title ?
	xlabel("Number of sub-elements");
	ylabel("Mass ratio [-]");
	grid;
    
	% Compute the relative differences between sdiv steps.
	reldiffs = abs(diff(freqs_serie)./freqs_serie(2:end, :));
    
    % Compute the relative differences between modes amplitudes
    mode_diff = zeros(numel(sdiv_serie)-1, BS.nbDOF, SOL.nbMode);
    norm_mode_diff = zeros(numel(sdiv_serie)-1, SOL.nbMode);

    for i = 1:numel(sdiv_serie)-1
        mode_diff(i, :, :) = abs(modes_serie(i+1, :, :) -  modes_serie(i, :, :));
        
        for m = 1:SOL.nbMode
            norm_mode_diff(i,m) = norm(mode_diff(i, :, m));
        end
    end
    
    % MAC
    % this section computes the MAC for each successive mode
    MAC_arr = zeros(numel(sdiv_serie)-1, SOL.nbMode, SOL.nbMode);
    
    for i = 1:numel(sdiv_serie)-1
        MAC_arr(i, :, :) = MAC(modes_serie(i, :, :), modes_serie(i+1, :, :), SOL);
    end
    
    % NOTE pretty certain MAC(mi, mi+1) should converge to eye(nbMode)
    % if modes do converge ?

    % Instantiate a figure object.
	figure("WindowStyle", "docked");
    
    % Plot of MAC convergence
    for i = 1:numel(sdiv_serie)-1
        subplot(1, numel(sdiv_serie), i);
        image(reshape(MAC_arr(i,:,:), SOL.nbMode, SOL.nbMode) ,'CDataMapping','scaled')
        colorbar
    end
    
    % Frequencies convergence
    figure("WindowStyle", "docked");
	% Plot frequencies w.r.t. sdiv.
	subplot(1, 3, 1);
	plot(sdiv_serie, freqs_serie);
	title("Frequencies convergence");
	xlabel("Number of sub-elements");
	ylabel("Natural frequency (Hz)");
	grid;

	% Plot relative differences w.r.t. sdiv.
	subplot(1, 3, 2);
	semilogy(sdiv_serie(2:end), reldiffs);
	title("Freqs reldiffs convergence");
	xlabel("Number of sub-elements");
	ylabel("Freqs relative difference");
	grid;
    
    % Mode convergence
    subplot(1, 3, 3);
	plot(sdiv_serie(2:end), norm_mode_diff);
	title("Mode L2 norm convergence");
	xlabel("Number of sub-elements");
	ylabel("Relative difference mode norm");
	grid;
end

 function mac = MAC(m1, m2, SOL)
    % MAC Compute Modal Assurance Criterion
    %
    % Arguments
    %   m1   1 x nbMode x nbDOF
    %   m2   1 x nbMode x nbDOF
    %
    % Returns
    %   mac  nbMode x nbMode
    
    % mac(m1, m2)[i][j] = ( (m1[i]' * m2[j]) / (|m1| * |m2|) )^2
    mac = zeros(SOL.nbMode, SOL.nbMode);
    for i = 1:SOL.nbMode
        for j = 1:SOL.nbMode

            norm_m1 = norm(m1(1, :, i));
            norm_m2 = norm(m2(1, :, j));

            mac(i,j) = (1 / (norm_m1 * norm_m2) * m1(1, :, i) * m2(1, :, j)').^2;
        end
    end
 end