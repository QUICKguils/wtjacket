function [frequencyHertz] = reduction(SdivStruct, AlgSys, nMode)
% TODO: documentation
% REDUCTION

ndList = SdivStruct.nodeList;
% Nodes highlighted in the structure
redNodesList = [18, 22];
redDOFsList = [];

for k = 1:length(redNodesList)
	redDOFsList = cat(2, redDOFsList, cat(2, ndList{1, redNodesList(k)}.dof(1:3), ndList{1, redNodesList(k)}.dof(6)));
end

condDOFsList = zeros(1, AlgSys.nDof-length(redDOFsList));
% TODO more elegant solution, this does work tho
k = 1;
l = 1;
while k <= AlgSys.nDof
	if not(ismember(k, redDOFsList))
		condDOFsList(l) = k;
		l = l + 1;
	end
	k = k + 1;
end

% Guyans-Irons Method - STATIC CONDENSATION
frequencyHertz = GuyanIronsReduction(AlgSys, nMode, redDOFsList, condDOFsList, redNodesList);

% Craig-Brampton Method
% TODO refactor this part into analysis/reduction_analysis
for m = 1:10
	frequencyHertz = [frequencyHertz CraigBramptonReduction(AlgSys, nMode, m, redDOFsList, condDOFsList, redNodesList)];
end
% TODO - freq convergence plot

% TODO - time integration
end

function [frequencyHertz] = GuyanIronsReduction(AlgSys, nMode, redDOFsList, condDOFsList, redNodesList)
% TODO refactor as suppressed collumns ?

% K submatrices.
KCC = AlgSys.K(condDOFsList, condDOFsList);
KCR = AlgSys.K(condDOFsList, redDOFsList);
KRR = AlgSys.K(redDOFsList, redDOFsList);
KRC = AlgSys.K(redDOFsList, condDOFsList);

% M Submatrices.
MCC = AlgSys.M(condDOFsList, condDOFsList);
MCR = AlgSys.M(condDOFsList, redDOFsList);
MRR = AlgSys.M(redDOFsList, redDOFsList);
MRC = AlgSys.M(redDOFsList, condDOFsList);

% Compute reduced AlgSys.
AlgSys.M = MRR - MRC / KCC * KCR - KRC / KCC * MCR + KRC / KCC * MCC / KCC * KCR;
AlgSys.K = KRR - KRC / KCC * KCR;

% Solve the eigenvalue problem.
[~, eigvals] = eigs(AlgSys.K, AlgSys.M, nMode, 'sm');
% Extract the natural frequencies.
frequencyRad  = sqrt(diag(eigvals));
frequencyHertz = frequencyRad / (2*pi);
end

function [frequencyHertz] = CraigBramptonReduction(AlgSys, nMode, m, redDOFsList, condDOFsList, redNodesList)

% K submatrices.
KII = AlgSys.K(condDOFsList, condDOFsList);
KIB = AlgSys.K(condDOFsList, redDOFsList);
KBB = AlgSys.K(redDOFsList, redDOFsList);
KBI = AlgSys.K(redDOFsList, condDOFsList);
% M submatrices.
MII = AlgSys.M(condDOFsList, condDOFsList);
MIB = AlgSys.M(condDOFsList, redDOFsList);
MBB = AlgSys.M(redDOFsList, redDOFsList);
MBI = AlgSys.M(redDOFsList, condDOFsList);

% Solve substructure fixed at the boundary - for m modes.
[XI, Omega] = eigs(KII, MII, m, 'sm');

% Mass-normalize modes per definition.
for i = 1:m
	XI(:, i) = XI(:, i) / sqrt(XI(:, i)' * MII * XI(:, i));
end

% Build reduced K.
KBBtilde = KBB - KBI / KII * KIB;
K_reduced = [KBBtilde zeros(length(redDOFsList), m); zeros(m, length(redDOFsList)) Omega];

% Build reduced M.
MBBtilde = MBB - (MBI/KII) * MIB + (KBI/KII) * (MII/KII) * KIB;
MBtilde = XI' * (MIB - (MII / KII) * KIB);
M_reduced = [MBBtilde MBtilde'; MBtilde eye(m)];

% Reduced AlgSys forr.
AlgSys.M = M_reduced;
AlgSys.K = K_reduced;

% Solve the reduced eigenvalue problem.
[~, eigvals] = eigs(AlgSys.K, AlgSys.M, nMode, 'sm');

% Extract the first natural frequencies.
frequencyRad = sqrt(diag(eigvals));
frequencyHertz = frequencyRad / (2*pi);
end
