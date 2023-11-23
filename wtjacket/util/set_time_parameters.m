%% Set the time discretization

function TimeParams = set_time_parameters(timeSample, initialConditions)
    % SET_TIME_PARAMETERS  Set the temporal parameters of the problem.
    
    % Ensure the time vector is in the expected shape.
    timeSample = reshape(timeSample, 1, []);
    
    % TODO: see if start and end fields are used.
    TimeParams.sample = timeSample;
    TimeParams.steps  = [diff(timeSample), timeSample(end)-timeSample(end-1)];
    TimeParams.start  = timeSample(1);
    TimeParams.end    = timeSample(end);
    TimeParams.numel  = numel(timeSample);
    TimeParams.initialConditions = initialConditions;
end