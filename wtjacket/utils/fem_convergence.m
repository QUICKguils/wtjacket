function fem_convergence
% CONVERGENCE  Analyze the convergence of the FEM simulation.

% TODO:
% - the solver (eigs) complain about ill conditioning from sdiv=4.
% - load the NX freqs with more precision
% - make an opts to write MCONV in MAT file ?
% - analyse the NX modes, I don't know how to load them.
% - implement other ways to analyze the convergence, e.g. MAC w.r.t. the
%   NX modes ?

MLSet = gen_matlab_sol();
NXSet = load_nx_sol();

plot_freq_convergence(MLSet);
plot_freq_convergence(NXSet);

plot_mode_convergence(MLSet);
plot_mode_convergence(NXSet);


end

function MLSet = gen_matlab_sol(sdiv_set)
% GEN_MATLAB_SOL  Generate the solutions for a set of subdivisions.
%
% Argument:
%	sdiv_set (1xN int) -- Set of desired number of subdivisions (default is 1:12).

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
if nargin == 0
	sdiv_set = 1:12;
end

% Gather the frequencies, modes, masses for the desired subdivisions.
freq_set = zeros(numel(sdiv_set), SOL.nbMode);
modes_set = zeros(numel(sdiv_set), BS.nbDOF, SOL.nbMode);    
masses_set = zeros(numel(sdiv_set), 2);

for i = 1:numel(sdiv_set)
	modeling(sdiv_set(i), '');
	SOL = load("res\modeling_sol.mat");
  BS = load(fullfile(res_file, "bare_struct.mat"));
  
  freq_set(i, :) = SOL.freqs;
  
  for d = 1:BS.nbDOF
      modes_set(i, d, :) = SOL.modes(d,:);
  end

  masses_set(i, 1) = SOL.mass_rbm;
  masses_set(i, 2) = SOL.mass_rbm/BS.mass;
  
end

% Build return data structure.
MLSet.sdiv_set = sdiv_set;
MLSet.freq_set = freq_set;
MLSet.modes_set = modes_set;
MLSet.masses_set = masses_set;
MLSet.name     = 'Matlab';
end

function NXSol = load_nx_sol
% LOAD_NX_SOL  Load the solution found via the NX simulations.

% Sequence of number of subdivisions.
NXSol.sdiv_set  = 1:8;

% WARN: 8-th mode for sdiv=1 is a 2nd order rotation.
NXSol.freq_set = [
	4.481593E-01, 4.553365E-01, 9.845646E-01, 6.926734E+00, 7.379443E+00, 1.656722E+01, 2.061273E+01, 2.249718E+01;
	4.481600E-01, 4.553359E-01, 9.845593E-01, 6.904420E+00, 7.353245E+00, 1.653022E+01, 1.993304E+01, 2.169000E+01;
	4.481601E-01, 4.553358E-01, 9.845590E-01, 6.902639E+00, 7.351166E+00, 1.652731E+01, 1.987232E+01, 2.161503E+01;
	4.481601E-01, 4.553358E-01, 9.845590E-01, 6.902199E+00, 7.350655E+00, 1.652662E+01, 1.985790E+01, 2.159727E+01;
	4.481601E-01, 4.553358E-01, 9.845589E-01, 6.902028E+00, 7.350457E+00, 1.652636E+01, 1.985248E+01, 2.159061E+01;
	4.481601E-01, 4.553358E-01, 9.845589E-01, 6.901943E+00, 7.350360E+00, 1.652624E+01, 1.984988E+01, 2.158741E+01;
	4.481601E-01, 4.553358E-01, 9.845589E-01, 6.901895E+00, 7.350305E+00, 1.652617E+01, 1.984842E+01, 2.158563E+01;
	4.481601E-01, 4.553358E-01, 9.845589E-01, 6.901866E+00, 7.350271E+00, 1.652612E+01, 1.984752E+01, 2.158453E+01];

NXSol.name = 'NX';
end

function plot_freq_convergence(SolSet)
% PLOT_FREQ_CONVERGENCE  Generate convergence graphs for the frequencies.
%
% Arguments:
%	SolSet (struct) -- Set of solution, with fields:
%	  sdiv_set (1xN int)    -- Set of number of subdivisions.
%	  freq_set (1xN double) -- Set of associated frequencies.
%	  name     (char)       -- Name of the solution.

% Build the solution name.
if SolSet.name ~= ""
	SolSet.name = " (" + SolSet.name + ")";
end

% Instantiate a figure object.
figure("WindowStyle", "docked");

% Plot the frequencies.
subplot(1, 2, 1);
plot(SolSet.sdiv_set, SolSet.freq_set);
title("Frequencies convergence" + SolSet.name);
xlabel("Number of sub-elements");
ylabel("Natural frequency (Hz)");
grid;

% Plot the residuals.
subplot(1, 2, 2);
% NOTE:
% this way of defining residuals is taken from the fluid mechanics
% course, section 3.3.
residuals = abs(diff(SolSet.freq_set)) ./ abs(SolSet.freq_set(2, :) - SolSet.freq_set(1, :));
semilogy(SolSet.sdiv_set(1:end-1), residuals);
title("Residuals convergence" + SolSet.name);
xlabel("Number of sub-elements");
ylabel("Residual");
grid;
end

function plot_mode_convergence(SolSet)
% TODO: implement this badboy
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

    % TODO mass convergence analysis

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
    
    % Mode convergence
    subplot(1, 3, 3);
	plot(sdiv_serie(2:end), norm_mode_diff);
	title("Mode L2 norm convergence");
	xlabel("Number of sub-elements");
	ylabel("Relative difference mode norm");
	grid;
end