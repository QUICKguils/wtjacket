function fem_convergence
	% CONVERGENCE  Analyze the convergence of the FEM simulation.

	% TODO:
	% - analyse the NX modes.
	% - implement other ways to analyze the convergence, e.g. MAC w.r.t. the
	%   NX modes ?

	% Directory where the present file lies.
	file_dir = fileparts(mfilename("fullpath"));
	% Results directory location.
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

	% Instantiate a figure object.
	figure("WindowStyle", "docked");

	% Plot frequencies.
	subplot(1, 2, 1);
	plot(sdiv_serie, freqs_serie);
	title("Frequencies convergence");
	xlabel("Number of sub-elements");
	ylabel("Natural frequency (Hz)");
	grid;

	% Plot relative differences w.r.t. the 8-th.
	subplot(1, 2, 2);
	reldiffs = abs(freqs_serie(1:end-1, :) - freqs_serie(end, :)) ./ freqs_serie(end, :);
	semilogy(sdiv_serie(1:end-1), reldiffs);
	title("Freqs reldiffs convergence");
	xlabel("Number of sub-elements");
	ylabel("Freqs relative difference");
	grid;

	% MAC matrix.

end

function NX = load_nx_sol
	% LOAD_NX_SOL Load the solution found via the NX simulation.

	NX.freq_1 = 4.481601E-01;
	NX.freq_2 = 4.553358E-01;
	NX.freq_3 = 9.845590E-01;
	NX.freq_4 = 6.902199E+00;
	NX.freq_5 = 7.350655E+00;
	NX.freq_6 = 1.652662E+01;
	NX.freq_7 = 1.985790E+01;

	NX.mode_1 = 0 ;
	NX.mode_2 = 0 ;
	NX.mode_3 = 0 ;
	NX.mode_4 = 0 ;
	NX.mode_5 = 0 ;
	NX.mode_6 = 0 ;
end
