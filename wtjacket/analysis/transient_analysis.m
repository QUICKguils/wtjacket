function transient_analysis()
% TRANSIENT_ANALYSIS  Analyze the results of the transient part.

% TODO:
% - stem plots of the phi(t)
% - frf(18, 18)
% - displ for the peak frequency of frf(18, 18).
% - frf of newmark's solution, both at the excitation point and at the
%   rotor location.

close all;

Cst = load_constants();
sdiv = 3;
nMode = 8;
[~, SdivStruct, AlgSys, FemSol] = modeling(Cst, sdiv, nMode, '');
[AlgSys, TransientSol] = transient(Cst, SdivStruct, AlgSys, FemSol, nMode, '');

analyze_phi(TransientSol);
% analyse_frf(AlgSys, FemSol, TransientSol.DiscreteLoad, SdivStruct.nodeList, TransientSol.ModalSup.mu);
analyse_NewmarkSol(SdivStruct, TransientSol)

end


function analyze_phi(TransientSol)
% plot of the phi's abs value.
%
% Analyse the modal participation factors.


phiAmplitude = max(TransientSol.ModalSup.phi, [], 2);

figure("WindowStyle", "docked");
stem(phiAmplitude);
set(gca,'yscal','log');
grid;
title([ ...
	"Distribution of modal participation factors (\phi)", ...
	"across the first " + num2str(TransientSol.ModalSup.nMode) + " modes"]);
xlabel("Associated mode number");
ylabel("max(\phi)");  % [m²/s² or 1/s²].

end

% TODO:
% - make this more flexible, not only (18, 18)
function analyse_frf(AlgSys, FemSol, DiscreteLoad, nodeList, mu)
% plot of the frf
% plot of the transient resp. via newmark for peak frequency.

	function frf = build_frf(AlgSys, FemSol, mu, w)
		frf = zeros(AlgSys.nDof_free);
		for s = 1:FemSol.nMode
			frf = frf + FemSol.mode(:, s) * FemSol.mode(:, s)' ./ (mu(s) * (FemSol.frequencyRad(s)^2 - w^2 + 2i*AlgSys.eps(s)*FemSol.frequencyRad(s)*w));
		end
	end

frf = @(w) build_frf(AlgSys, FemSol, mu, w);

freq = (0:0.01:20)*(2*pi);
nFreq = numel(freq);
xSample = zeros(AlgSys.nDof_free, nFreq);

for iFreq = 1:nFreq
	xSample(:, iFreq) = abs(frf(freq(iFreq)) * DiscreteLoad.spatial);
end

x_proj = xSample(nodeList{18}.dof(1), :) * cosd(45) + xSample(nodeList{18}.dof(2), :) * sind(45);

figure("WindowStyle", "docked");
semilogy(freq/(2*pi), x_proj);
grid;
xlabel("Frequency (Hz)");
ylabel("Displacement (node 18, dir of load)");
title("Amplitude of the steady harmonic load");
% -> pic vers 4e et 5e comme attendu. Mais en fait plus gros pic encore
% vers 0.44Hz, càd vers la première fréq. Sion essaye de changer la
% fréquence d'excitation de 1hz à celle-ci, on voit qu'en effet la
% réponse transitoire s'amplifie fort au cours du temps, jusqu'à
% atteindre une valeur fort grande.
end


% WARN: the fft strongly depend on the chosen time interval.
% make sure to integrate enough 1st mode periods.
function analyse_NewmarkSol(SdivStruct, TransientSol)
% Plot of the FFT of the Newmark solution.

y = fft(TransientSol.Newmark.q(SdivStruct.nodeList{18}.dof(1), :));
fs = 1/TransientSol.TimeParams.steps(1);
% f = (0:length(y)-1)*fs/length(y);
n = length(TransientSol.Newmark.q(SdivStruct.nodeList{18}.dof(1), :));
fshift = (-n/2:n/2-1)*(fs/n);
yshift = fftshift(y);
figure("WindowStyle", "docked");
semilogy(fshift, abs(yshift)); grid;
semilogy(fshift(ceil(end/2:end)), abs(yshift(ceil(end/2:end))));
title("FFT of the transient response");
xlabel("Frequency (Hz)");
ylabel("FFT");
xlim([0, 30]);
grid;

end
