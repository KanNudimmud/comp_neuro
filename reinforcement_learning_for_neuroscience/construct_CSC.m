function x = construct_CSC(s)
    % Create complete serial compound representation.
    %
    % USAGE: x = construct_CSC(s)
    %
    % INPUTS:
    %   s - stimulus timeseries (see construct_stimulus.m)
    %
    % OUTPUTS:
    %   x - [trial_length x trial_length*D] complete serial compound representation
    %
    % Credit: Sam Gershman
    
    x = [];
    for d = 1:size(s,2)
        x = [x diag(s(:,d))];
    end