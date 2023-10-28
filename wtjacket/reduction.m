function varargout = reduction(opts)
% REDUCTION  

% X. Gather and return the relevant calculated data

optrets = {AlgSys, TransientSol};
varargout(1:nargout) = optrets(1:nargout);

end