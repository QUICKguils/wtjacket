function convergence
	% CONVERGENCE  Analyze the convergence of the modeling solution.

	% Directory where the present file lies.
	file_dir = fileparts(mfilename("fullpath"));
	% Results file.
	res_file = (fullfile(file_dir, "../../res/"));
	
	% Run the modeling once, if SOL data are not available.
	if ~isfile(fullfile(res_file, "modeling_sol.mat"))
		% Perform a FEM simulation.
		%
		% 1. Reset class internal states, close previous plots.
		clear Node Elem;
		close all;
		% 2. Initialize MAT file.
		constants();
		bare_struct('');
		% 3. Run the FEM simulation.
		modeling(1, '');
	end
	SOL = load(fullfile(res_file, "modeling_sol.mat"));

	% Desired sdiv sequence.
	sdiv_serie = 1:8;
	freqs_serie = zeros(numel(sdiv_serie), SOL.nbMode);

	% Gather the frequencies for the desired subdivisions.
	for i = 1:numel(sdiv_serie)
		% Perform a FEM simulation.
		%
		% 1. Reset class internal states, close previous plots.
		clear Node Elem;
		close all;
		% 2. Initialize MAT file.
		bare_struct('');
		% 3. Run the FEM simulation.
		modeling(sdiv_serie(i), '');
		SOL = load("res\modeling_sol.mat");
		freqs_serie(i, :) = SOL.freqs;
	end

	% Compute the relative differences bw. sdiv steps.
	reldiffs = abs(diff(freqs_serie)./freqs_serie(2:end, :));

	% Plot this badboy bruh.
	semilogy(sdiv_serie(2:end), reldiffs);
	title("reldiffs");
	xlabel("Number of subelements");
	ylabel("Relative difference");
	grid;

	% Plot the freqs.
	figure;
	plot(sdiv_serie, freqs_serie);
	title("Freqs");
	xlabel("Number of subelements");
	ylabel("Freqs (Hz)");
	grid;
end
