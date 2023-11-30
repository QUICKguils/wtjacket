function reduction_analysis(varargin)
% REDUCTION_ANALYSIS  Analyze the results of the reduction part.
%
% Argument:
%	mSet (1xN int) -- Number of modes in the reduction (default: 1:50).

%% Set input arguments

% Set default value for optional inputs.
optargs = {1:50};
% Overwrite default value of optional inputs.
optargs(1:numel(varargin)) = varargin;
% Place optional args in memorable variable names.
[mSet] = optargs{:};

%% Precompute dependencies

% Fetch and overwrite the defaults execution parameters.
RunArg = load_defaults();
RunArg.opts = '';
% Fetch the project statement data.
Stm = load_statement();
% Compute the FEM solution.
[~, SdivStruct, AlgSys, FemSol] = modeling(RunArg, Stm);
% Compute the mode superposition parameters.
[AlgSys, ~] = transient(RunArg, Stm, SdivStruct, AlgSys, FemSol);


%% Compute convergence data

% Taking ref frequencies*1.02 as target.
ref_freq = FemSol.frequencyHertz;
target_ref_freq = 1.02 .* ref_freq;

% For m = 0, kinda.
tic;
RunArg.m = 1;
ReductionSol = reduction(RunArg, Stm, SdivStruct, AlgSys);
GIReducedFemSol = ReductionSol.GIReducedFemSol;
t_GI = toc;
freq_array = [GIReducedFemSol.frequencyHertz];


% Non-exhaustive reduced system study.
n_m = length(mSet);
t_CB = zeros(1, n_m);

for i=1:n_m
	tic;
	RunArg.m = mSet(i);
	ReductionSol = reduction(RunArg, Stm, SdivStruct, AlgSys);
	CBReducedFemSol = ReductionSol.CBReducedFemSol;
	t_CB(i) = toc;
	freq_array = [freq_array CBReducedFemSol.frequencyHertz];
end

mSet = [0 mSet];
t_CB = [t_GI t_CB];
t_CB = 1000 * t_CB;

% relative error
rel_error_f = 100 * abs((freq_array-target_ref_freq)./target_ref_freq);

%% Plot the convergence results

% plot freq convergence
figure("WindowStyle", "docked");
plot(mSet, freq_array);
xlabel("Number of sub-eigenmodes taken for R matrix computation");
ylabel("Frequency (Hz)");
ylim([0, 100])
grid;

% rel convergence
figure("WindowStyle", "docked");
plot(mSet, rel_error_f);
xlabel("Number of sub-eigenmodes taken for R matrix computation");
ylabel("% Absolute Relative Error on Frequencies (-)");
ylim([0, 10])
grid;

% plot computation time
figure("WindowStyle", "docked");
plot(mSet, t_CB);
xlabel("Number of sub-eigenmodes taken for R matrix computation");
ylabel("Computation time (ms)");
grid;

% plot displacements : approx wrt. exact 
figure("WindowStyle", "docked");

ThisLoad = SdivStruct.loadList{1};
loadDirection = ThisLoad.direction;

ReductionSol = reduction(RunArg, Stm, SdivStruct, AlgSys);
ReducedNewmarkSol = ReductionSol.ReducedNewmarkSol;
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
