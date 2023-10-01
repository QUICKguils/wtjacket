function fem_convergence
	% CONVERGENCE  Analyze the convergence of the FEM simulation.

	% TODO:
	% - analyse the NX modes.
	% - implement other ways to analyze the convergence, e.g. MAC w.r.t. the
	%   NX modes ?

	% Directory where the present file lies.
	file_dir = fileparts(mfilename("fullpath"));
	% Results file.
	res_file = (fullfile(file_dir, "../../res/"));

	% Run the bare structure simulation once, if SOL data are not available.
	if ~isfile(fullfile(res_file, "modeling_sol.mat"))
		modeling(1, '');
	end
	SOL = load(fullfile(res_file, "modeling_sol.mat"));

	% Desired sdiv sequence.
	sdiv_serie = 1:8;
	freqs_serie = zeros(numel(sdiv_serie), SOL.nbMode);

	% Gather the frequencies for the desired subdivisions.
	for i = 1:numel(sdiv_serie)
		modeling(sdiv_serie(i), '');
		SOL = load("res\modeling_sol.mat");
		freqs_serie(i, :) = SOL.freqs;
	end

	% Compute the relative differences between sdiv steps.
	reldiffs = abs(diff(freqs_serie)./freqs_serie(2:end, :));

	% Instantiate a figure object.
	figure("WindowStyle", "docked");

	% Plot frequencies w.r.t. sdiv.
	subplot(1, 2, 1);
	plot(sdiv_serie, freqs_serie);
	title("Frequencies convergence");
	xlabel("Number of sub-elements");
	ylabel("Natural frequency (Hz)");
	grid;

	% Plot relative differences w.r.t. sdiv.
	subplot(1, 2, 2);
	semilogy(sdiv_serie(2:end), reldiffs);
	title("Freqs reldiffs convergence");
	xlabel("Number of sub-elements");
	ylabel("Freqs relative difference");
	grid;
end
