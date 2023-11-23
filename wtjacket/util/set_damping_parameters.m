%% Set the proportional damping parameters

function [C, eps] = set_damping_parameters(eps, w0, AlgSys)
    % SET_DAMPING_MATRIX  Set the damping matrix, assuming a proportional damping.
    %
    % Arguments:
    %	eps (1x2 double)
    %	  Proportional damping ratios of the first two modes.
    %	w0 (1xnMode double)
    %	  Natural frequencies of the asociated conservative system [rad/s].
    %	AlgSys (struct)
    %	  Parameters of the discrete algebraic system.
    % Returns:
    %	C   (NxN double)     -- Proportional damping matrix.
    %	eps (1xNmode double) -- Proportional damping ratios of the first modes.
    %
    % See reference book, p.156.
    
    a = 2             * (w0(1)*eps(1) - w0(2)*eps(2)) / (w0(1)^2 - w0(2)^2);
    b = 2*w0(1)*w0(2) * (w0(1)*eps(2) - w0(2)*eps(1)) / (w0(1)^2 - w0(2)^2);
    
    C   = a*AlgSys.K + b*AlgSys.M;
    eps = 0.5 * (a*w0 + b./w0);
end