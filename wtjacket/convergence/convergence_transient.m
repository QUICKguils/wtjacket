function convergence_transient(TransientSol, SdivStruct)
% CONVERGENCE_TRANSIENT  Analyze the convergence of the transient solutions.

% TODO:
% Some of these plots don't analyze the convergence of the solution, but
% just the nature of the solution. These plots should be put in a
% dedicated section, in the main transient.m file.

% Analyse the modal participation factors.
% We see that 
figure("WindowStyle", "docked");
phiAmplitude = max(TransientSol.ModalSup.phi, [], 2);
stem(phiAmplitude);
set(gca,'yscal','log');
% NOTE: not well expressive on a ylog axis.
grid;
title("Modal participation factors");

y = fft(TransientSol.Newmark.q(SdivStruct.nodeList{18}.dof(1), :));
fs = 1/TransientSol.TimeParams.steps(1);
f = (0:length(y)-1)*fs/length(y);
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