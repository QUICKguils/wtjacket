function reduction_analysis(Stm, SdivStruct, AlgSys, FemSol)
% REDUCTION_ANALYSIS  Analyze the results of the reduction part.

% Packing relevant execution parameters.
RunArg.sdiv = 0;
RunArg.nMode = 8;
RunArg.opts = '';

% taking ref frequencies*1.02 as target
ref_freq = FemSol.frequencyHertz;
target_ref_freq = 1.02 .* ref_freq;

% for m = 0, kinda
tic;
[~, ~, GIReducedFemSol, ~, ~, ~] = reduction(RunArg, Stm, SdivStruct, AlgSys, 1);
t_GI = toc
freq_array = [GIReducedFemSol.frequencyHertz];


% non-exhaustive reduced system study
m_arr = 1:50;
n_m = length(m_arr)
t_CB = zeros(1, n_m);

for i=1:n_m
	tic;
	[~, ~, ~, ~, CBReducedFemSol, ~] = reduction(RunArg, Stm, SdivStruct, AlgSys, m_arr(i));
	t_CB(i) = toc;
	freq_array = [freq_array CBReducedFemSol.frequencyHertz];
end

m_arr = [0 m_arr];
t_CB = [t_GI t_CB];
t_CB = 1000 * t_CB;

% relative error 
rel_error_f = 100 * abs((freq_array-target_ref_freq)./target_ref_freq);

% plot freq convergence
figure("WindowStyle", "docked");
plot(m_arr, freq_array);
xlabel("Number of sub-eigenmodes taken for R matrix computation");
ylabel("Frequency (Hz)");
ylim([0, 100])
grid;

% rel convergence
figure("WindowStyle", "docked");
plot(m_arr, rel_error_f);
xlabel("Number of sub-eigenmodes taken for R matrix computation");
ylabel("% Absolute Relative Error on Frequencies (-)");
ylim([0, 10])
grid;

% plot computation time
figure("WindowStyle", "docked");
plot(m_arr, t_CB);
xlabel("Number of sub-eigenmodes taken for R matrix computation");
ylabel("Computation time (ms)");
grid;

end