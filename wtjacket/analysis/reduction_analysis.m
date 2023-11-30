function reduction_analysis(Stm, SdivStruct, AlgSys, FemSol, NewmarkSol, ReducedNewmarkSol)
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
m_arr = 1:25;
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

% plot displacements : approx wrt. exact 
figure("WindowStyle", "docked");

ThisLoad = SdivStruct.loadList{1};
loadDirection = ThisLoad.direction;

[GIReducedSdivStruct, GIReducedAlgSys, GIReducedFemSol, CBReducedAlgSys, CBReducedFemSol, ReducedNewmarkSol] = reduction(RunArg, Stm, SdivStruct, AlgSys, 3);
t = ReducedNewmarkSol.TimeParams.sample;

% upper graph
q = NewmarkSol.q;
qX = q(103, :) * loadDirection(1);
qY = q(104, :) * loadDirection(2);
qZ = q(105, :) * loadDirection(3);
qDir = qX + qY + qZ;

qr = ReducedNewmarkSol.q;
qrX = qr(1, :) * loadDirection(1);
qrY = qr(2, :) * loadDirection(2);
qrZ = qr(3, :) * loadDirection(3);
qrDir = qrX + qrY + qrZ;

subplot(2, 1, 1);
plot(t, 100*(qDir - qrDir)./qDir);
grid;

% lower graph
q = NewmarkSol.q;
qX = q(127, :) * loadDirection(1);
qY = q(128, :) * loadDirection(2);
qZ = q(129, :) * loadDirection(3);
qDir = qX + qY + qZ;

qr = ReducedNewmarkSol.q;
qrX = qr(5, :) * loadDirection(1);
qrY = qr(6, :) * loadDirection(2);
qrZ = qr(7, :) * loadDirection(3);
qrDir = qrX + qrY + qrZ;

subplot(2, 1, 2);
plot(t, 100*(qDir - qrDir)./qDir);
grid;
end